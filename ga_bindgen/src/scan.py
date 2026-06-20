"""Scan ga/*.hpp with libclang and emit a manifest.json.

Usage:  python3 ga_bindgen/src/scan.py [--out PATH]

Produces a structured manifest matching `model.Manifest`. Symbols outside
`ga/` (system headers, fmt, etc.) are filtered out.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from collections import defaultdict
from pathlib import Path

import clang.cindex as cx

# Allow running as `python ga_bindgen/src/scan.py` from anywhere.
sys.path.insert(0, str(Path(__file__).parent))
from clang_setup import GA_DIR, PROJECT_ROOT, compile_args, in_ga_tree

_OUTER_NS_RE = re.compile(r"^hd::ga::(?=[A-Z_])")


def _rel_source(file) -> str:
    """Return source file path relative to the project root as a POSIX string.

    libclang returns absolute OS-native paths which differ between platforms.
    Storing a project-root-relative POSIX path makes manifest.json identical
    across macOS, Linux, and Windows regardless of where the repo is cloned.
    """
    if file is None:
        return ""
    try:
        return Path(file.name).relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return Path(file.name).as_posix()  # outside project root — keep absolute


def _norm_canonical(spelling: str) -> str:
    """Normalise canonical_underlying to be platform-independent.

    libclang on macOS spells template-instantiation underlying types with a
    fully-qualified outer name (e.g. 'hd::ga::Vec3_t<...>'), while libclang
    on Windows omits the leading namespace for such types ('Vec3_t<...>').
    Named class types in sub-namespaces (e.g. 'hd::ga::pga::Point2d<double>')
    are consistent on both platforms and are left unchanged.
    Strip the leading 'hd::ga::' only when it immediately precedes an
    uppercase letter — i.e. a PascalCase template name, not a sub-namespace.
    """
    return _OUTER_NS_RE.sub("", spelling)


from model import (
    SCHEMA_VERSION,
    Constant,
    Constructor,
    Enum,
    Field,
    FunctionGroup,
    Manifest,
    Overload,
    TypeAlias,
)

# Headers we parse to build the unified manifest. ga_pga.hpp pulls in the
# pga side; ga_ega.hpp pulls in the ega side. ga_sta.hpp pulls in the sta side.
# We parse all and merge.
SOURCE_HEADERS = ["ga/ga_ega.hpp", "ga/ga_pga.hpp", "ga/ga_sta.hpp"]

# Namespaces we emit bindings for. `hd::ga::detail` is internal — excluded.
TARGET_NAMESPACES = {"hd::ga", "hd::ga::ega", "hd::ga::pga", "hd::ga::sta"}

# Namespace-scope types the structural rules below would otherwise pick up, but which are
# INTERNALS of the stateful (unbound) dynamic_system{2,3}dp force / contact / integration
# machinery. The systems themselves are deliberately not bound — Python reconstructs their
# behaviour from the bound primitives (see ga_py/README.md) — so these helper structs, the
# std::function callback aliases, and the integrator selector are not independently usable
# from Python. Skipping them by name keeps the binding aligned with that design and prevents a
# future regeneration from surprise-binding them (or choking on the std::function aliases,
# which are not GA value types). Unlike the bound PODs (pose/kin_state/joint_state), none of
# these is user-facing data a Python user would construct and pass around.
SKIP_TYPE_SPELLINGS = {
    "grounded_spring2dp",  # add_grounded_spring component bundle (dynamic_system internal)
    "grounded_spring3dp",
    "contact_state3dp",    # live state handed to a grinding force law (3D contact loop)
    "grinding_law3dp",     # std::function<vec3dp(contact_state3dp const&)> -- force-law callback
    "wrench_fn",           # std::function<bivecNdp(value_t)> -- applied-wrench callback
    "integrator_kind",     # selector for dynamic_system3dp::set_integrator (rk4 | abm2)
}


def _cursor_kind(cursor: cx.Cursor):
    """Return cursor.kind, or None if its kind id is unknown to these bindings.

    The pip `libclang` package lags LLVM and ships an incomplete CursorKind
    table, so a C++ construct newer than that table (e.g. the C++20
    parenthesized aggregate-init expression, kind id 155, emitted for member
    initialisers like `pose(origin, angle)`) makes CursorKind.from_id raise
    ValueError. Such cursors only ever occur inside function / initializer
    bodies, which the scanner never collects, so treating an unknown kind as
    "skip this cursor and its subtree" is safe and keeps the scan working
    against any libclang version.
    """
    try:
        return cursor.kind
    except ValueError:
        return None


def fq_namespace(cursor: cx.Cursor) -> str:
    parts: list[str] = []
    p = cursor.semantic_parent
    while p is not None and p.kind != cx.CursorKind.TRANSLATION_UNIT:
        if p.kind == cx.CursorKind.NAMESPACE:
            parts.append(p.spelling)
        p = p.semantic_parent
    return "::".join(reversed(parts))


def build_class_template_index(tu: cx.TranslationUnit) -> dict[str, cx.Cursor]:
    """Index every CLASS_TEMPLATE in ga/ by its short name (e.g. "Vec3_t")."""
    idx: dict[str, cx.Cursor] = {}

    def walk(c: cx.Cursor):
        k = _cursor_kind(c)
        if k is None:
            return
        if c.location.file and not in_ga_tree(c.location):
            if k != cx.CursorKind.NAMESPACE:
                return
        if k == cx.CursorKind.CLASS_TEMPLATE and c.spelling:
            idx.setdefault(c.spelling, c)
        for child in c.get_children():
            walk(child)

    walk(tu.cursor)
    return idx


_TAG_RE = re.compile(r"<\s*[^,<>]+?\s*,\s*([^,<>]+?)\s*>\s*$")


def extract_tag(type_spelling: str) -> str | None:
    """Pull the tag from a 2-arg template spelling.

    >>> extract_tag("hd::ga::MVec8_t<T, mvec3d_tag>")
    'mvec3d_tag'
    >>> extract_tag("hd::ga::MVec8_t<double, hd::ga::mvec3dp_e_tag>")
    'mvec3dp_e_tag'
    """
    m = _TAG_RE.search(type_spelling)
    if not m:
        return None
    return m.group(1).strip().rsplit("::", 1)[-1]


def build_partial_spec_index(
    tu: cx.TranslationUnit,
) -> dict[tuple[str, str], cx.Cursor]:
    """Index every CLASS_TEMPLATE_PARTIAL_SPECIALIZATION by (primary_name, tag)."""
    idx: dict[tuple[str, str], cx.Cursor] = {}

    def walk(c: cx.Cursor):
        k = _cursor_kind(c)
        if k is None:
            return
        if c.location.file and not in_ga_tree(c.location):
            if k != cx.CursorKind.NAMESPACE:
                return
        if k == cx.CursorKind.CLASS_TEMPLATE_PARTIAL_SPECIALIZATION:
            tag = extract_tag(c.type.spelling)
            if c.spelling and tag:
                idx.setdefault((c.spelling, tag), c)
        for child in c.get_children():
            walk(child)

    walk(tu.cursor)
    return idx


def template_name(canonical_spelling: str) -> str:
    """Extract the unqualified template name from a canonical type spelling.

    >>> template_name("hd::ga::Vec3_t<double, hd::ga::vec3d_tag>")
    'Vec3_t'
    """
    head = canonical_spelling.split("<", 1)[0]  # "hd::ga::Vec3_t"
    return head.rsplit("::", 1)[-1]  # "Vec3_t"


def _collect_ctors_and_fields(
    template_or_spec: cx.Cursor, fields_out: list[Field], ctors_out: list[Constructor]
) -> None:
    """Append public field decls and constructor signatures from a class
    template / partial specialization cursor."""
    for child in template_or_spec.get_children():
        if child.access_specifier != cx.AccessSpecifier.PUBLIC:
            continue
        if child.kind == cx.CursorKind.FIELD_DECL:
            fields_out.append(Field(name=child.spelling, type=child.type.spelling))
        elif child.kind == cx.CursorKind.CONSTRUCTOR:
            params = [
                p.type.spelling
                for p in child.get_children()
                if p.kind == cx.CursorKind.PARM_DECL
            ]
            ctors_out.append(Constructor(param_types=params))


def _detect_base_class(template_cursor: cx.Cursor) -> str | None:
    """If the class template has a single public base, return its short
    template name (e.g. 'Vec2d'). Used for PGA primitives like
    `Point2d : public Vec2d<T>`."""
    for child in template_cursor.get_children():
        if child.kind == cx.CursorKind.CXX_BASE_SPECIFIER:
            spelling = child.type.spelling  # e.g. "Vec2d<T>"
            if "<" in spelling:
                spelling = spelling.split("<", 1)[0]
            return spelling.strip().rsplit("::", 1)[-1]
    return None


def extract_type(
    typedef_cursor: cx.Cursor,
    template_index: dict[str, cx.Cursor],
    partial_spec_index: dict[tuple[str, str], cx.Cursor],
) -> TypeAlias:
    ns = fq_namespace(typedef_cursor)
    underlying = typedef_cursor.underlying_typedef_type.spelling
    canonical = typedef_cursor.underlying_typedef_type.get_canonical().spelling

    fields: list[Field] = []
    ctors: list[Constructor] = []
    base_class: str | None = None

    primary_name = template_name(canonical)
    template = template_index.get(primary_name)
    if template is not None and template.kind == cx.CursorKind.CLASS_TEMPLATE:
        _collect_ctors_and_fields(template, fields, ctors)
        base_class = _detect_base_class(template)

    # The user typedef may resolve to a partial specialization that adds
    # grade-composition constructors on top of the primary template's ctors
    # (e.g. MVec8_t<T, mvec3d_tag> adds MVec8_t(Scalar3d<T>, Vec3d<T>, ...)).
    tag = extract_tag(canonical)
    if tag:
        spec = partial_spec_index.get((primary_name, tag))
        if spec is not None:
            # Specialization may also re-declare fields, but for the shapes
            # we currently support, the fields come from the primary template.
            spec_ctors: list[Constructor] = []
            _collect_ctors_and_fields(spec, [], spec_ctors)
            ctors.extend(spec_ctors)

    src = typedef_cursor.location.file
    return TypeAlias(
        namespace=ns,
        name=typedef_cursor.spelling,
        underlying=underlying,
        canonical_underlying=_norm_canonical(canonical),
        fields=fields,
        constructors=ctors,
        base_class=base_class,
        source_file=_rel_source(src),
        source_line=typedef_cursor.location.line,
    )


def _is_pure_data_struct(cursor: cx.Cursor) -> bool:
    """True for an aggregate of only public data members — no methods, bases,
    constructors, or non-public members.

    This is exactly the shape of the physics POD structs (pose2dp/3dp,
    kin_state2dp/3dp) and deliberately excludes the stateful system classes
    (static_frame2dp, static_system2dp, kinematic_system2dp, …), which carry
    methods, inheritance, and private std::vector / std::string members.
    """
    children = list(cursor.get_children())
    if not children:
        return False
    for ch in children:
        if _cursor_kind(ch) != cx.CursorKind.FIELD_DECL:
            return False
        if ch.access_specifier != cx.AccessSpecifier.PUBLIC:
            return False
    return True


def extract_data_struct(cursor: cx.Cursor) -> TypeAlias:
    """Build a TypeAlias for a pure-data struct (see _is_pure_data_struct).

    Unlike the template-instantiation typedefs, these are concrete non-template
    structs, so `underlying` / `canonical_underlying` are just the qualified
    struct name and `constructors` is left empty — the emitter synthesises an
    aggregate constructor from the fields.
    """
    ns = fq_namespace(cursor)
    fields = [
        Field(name=ch.spelling, type=ch.type.spelling)
        for ch in cursor.get_children()
        if _cursor_kind(ch) == cx.CursorKind.FIELD_DECL
    ]
    src = cursor.location.file
    return TypeAlias(
        namespace=ns,
        name=cursor.spelling,
        underlying=cursor.spelling,
        canonical_underlying=f"{ns}::{cursor.spelling}",
        fields=fields,
        constructors=[],
        base_class=None,
        source_file=_rel_source(src),
        source_line=cursor.location.line,
    )


def extract_function(cursor: cx.Cursor) -> tuple[str, str, Overload]:
    """Return (namespace, short_name, Overload).

    Uses get_children()/PARM_DECL because cursor.get_arguments() returns
    empty for FUNCTION_TEMPLATE cursors in this libclang Python binding.
    """
    ns = fq_namespace(cursor)
    name = cursor.spelling
    return_type = cursor.result_type.spelling
    param_types = [
        p.type.spelling
        for p in cursor.get_children()
        if p.kind == cx.CursorKind.PARM_DECL
    ]
    src = cursor.location.file
    return (
        ns,
        name,
        Overload(
            return_type=return_type,
            param_types=param_types,
            source_file=_rel_source(src),
            source_line=cursor.location.line,
        ),
    )


def parse_header(header_path: Path, args: list[str]):
    idx = cx.Index.create()
    tu = idx.parse(
        str(header_path),
        args=args,
        options=cx.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD,
    )
    errors = [d for d in tu.diagnostics if d.severity >= cx.Diagnostic.Error]
    return tu, errors


def collect(
    tu: cx.TranslationUnit,
) -> tuple[
    list[TypeAlias],
    dict[tuple[str, str], list[Overload]],
    dict[tuple[str, str], list[Overload]],
    set[str],
    list[Constant],
    list[Enum],
]:
    """Walk the TU and classify cursors into types / functions / operators /
    namespaces / constants / enums."""
    types: dict[tuple[str, str], TypeAlias] = {}
    functions: dict[tuple[str, str], list[Overload]] = defaultdict(list)
    operators: dict[tuple[str, str], list[Overload]] = defaultdict(list)
    namespaces: set[str] = set()
    constants_index: dict[tuple[str, str], Constant] = {}
    enums_index: dict[tuple[str, str], Enum] = {}

    template_index = build_class_template_index(tu)
    partial_spec_index = build_partial_spec_index(tu)

    seen_overload_keys: set[tuple] = set()  # dedup across multiple parses

    def walk(cursor: cx.Cursor):
        k = _cursor_kind(cursor)
        if k is None:
            return
        if cursor.location.file and not in_ga_tree(cursor.location):
            # Don't recurse into headers outside ga/ — except for top-level
            # NAMESPACE cursors, since the same namespace spans multiple files.
            if k != cx.CursorKind.NAMESPACE:
                return

        if k == cx.CursorKind.NAMESPACE:
            ns = fq_namespace(cursor)
            full = f"{ns}::{cursor.spelling}" if ns else cursor.spelling
            namespaces.add(full)
        elif k == cx.CursorKind.TYPE_ALIAS_DECL and in_ga_tree(cursor.location):
            ns = fq_namespace(cursor)
            if ns in TARGET_NAMESPACES and cursor.spelling not in SKIP_TYPE_SPELLINGS:
                key = (ns, cursor.spelling)
                if key not in types:
                    types[key] = extract_type(
                        cursor, template_index, partial_spec_index
                    )
        elif (
            k in (cx.CursorKind.STRUCT_DECL, cx.CursorKind.CLASS_DECL)
            and in_ga_tree(cursor.location)
            and cursor.is_definition()
        ):
            # Concrete pure-data structs (pose2dp/3dp, kin_state2dp/3dp). Only
            # namespace-scope aggregates are bound; nested structs (e.g. inside
            # kinematic_system2dp) and the stateful system classes are skipped.
            parent = cursor.semantic_parent
            if (
                parent is not None
                and parent.kind == cx.CursorKind.NAMESPACE
                and fq_namespace(cursor) in TARGET_NAMESPACES
                and cursor.spelling not in SKIP_TYPE_SPELLINGS
                and _is_pure_data_struct(cursor)
            ):
                ns = fq_namespace(cursor)
                key = (ns, cursor.spelling)
                if key not in types:
                    types[key] = extract_data_struct(cursor)
        elif (
            k == cx.CursorKind.ENUM_DECL
            and in_ga_tree(cursor.location)
            and cursor.is_definition()
        ):
            # Namespace-scope (scoped or plain) enums in the target namespaces, e.g.
            # constraint2dp / joint2dp. Bound as nb::enum_; also unlocks pure-data
            # structs that carry an enum-typed field (loop_constraint2dp,
            # joint_state2dp). Nested enums (semantic parent not a NAMESPACE) are
            # skipped, mirroring the struct rule above.
            parent = cursor.semantic_parent
            if (
                parent is not None
                and parent.kind == cx.CursorKind.NAMESPACE
                and fq_namespace(cursor) in TARGET_NAMESPACES
                and cursor.spelling not in SKIP_TYPE_SPELLINGS
            ):
                ns = fq_namespace(cursor)
                key = (ns, cursor.spelling)
                if key not in enums_index:
                    enumerators = [
                        ch.spelling
                        for ch in cursor.get_children()
                        if _cursor_kind(ch) == cx.CursorKind.ENUM_CONSTANT_DECL
                    ]
                    src = cursor.location.file
                    # is_scoped_enum() may be absent on older libclang bindings;
                    # it is informational only (emission uses `Type::value` for both
                    # scoped and plain enums in C++11+), so default to True.
                    is_scoped = getattr(cursor, "is_scoped_enum", lambda: True)()
                    enums_index[key] = Enum(
                        namespace=ns,
                        name=cursor.spelling,
                        canonical=f"{ns}::{cursor.spelling}",
                        enumerators=enumerators,
                        scoped=is_scoped,
                        source_file=_rel_source(src),
                        source_line=cursor.location.line,
                    )
        elif k == cx.CursorKind.VAR_DECL and in_ga_tree(cursor.location):
            # Only namespace-scope constants — not locals inside functions.
            # libclang reports VAR_DECL for both, but the semantic parent
            # of a true namespace-scope constant is a NAMESPACE cursor.
            parent = cursor.semantic_parent
            if parent is not None and parent.kind == cx.CursorKind.NAMESPACE:
                ns = fq_namespace(cursor)
                if ns in TARGET_NAMESPACES:
                    key = (ns, cursor.spelling)
                    if key not in constants_index:
                        src = cursor.location.file
                        constants_index[key] = Constant(
                            namespace=ns,
                            name=cursor.spelling,
                            type=cursor.type.spelling,
                            source_file=_rel_source(src),
                            source_line=cursor.location.line,
                        )
        elif k in (cx.CursorKind.FUNCTION_DECL, cx.CursorKind.FUNCTION_TEMPLATE):
            if not in_ga_tree(cursor.location):
                pass
            else:
                ns, name, ov = extract_function(cursor)
                if ns in TARGET_NAMESPACES:
                    # Dedup: same overload may appear when both headers parsed.
                    dedup_key = (
                        ns,
                        name,
                        ov.return_type,
                        tuple(ov.param_types),
                        ov.source_file,
                        ov.source_line,
                    )
                    if dedup_key in seen_overload_keys:
                        return
                    seen_overload_keys.add(dedup_key)

                    if name.startswith("operator"):
                        operators[(ns, name)].append(ov)
                    else:
                        functions[(ns, name)].append(ov)

        for child in cursor.get_children():
            walk(child)

    walk(tu.cursor)
    return (
        list(types.values()),
        functions,
        operators,
        namespaces,
        list(constants_index.values()),
        list(enums_index.values()),
    )


def to_groups(d: dict[tuple[str, str], list[Overload]]) -> list[FunctionGroup]:
    return sorted(
        [
            FunctionGroup(namespace=ns, name=name, overloads=ovs)
            for (ns, name), ovs in d.items()
        ],
        key=lambda g: (g.namespace, g.name),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out",
        default=str(PROJECT_ROOT / "ga_bindgen" / "manifest.json"),
        help="Path to write manifest.json",
    )
    args_ns = ap.parse_args()

    args = compile_args()

    all_types: dict[tuple[str, str], TypeAlias] = {}
    all_funcs: dict[tuple[str, str], list[Overload]] = defaultdict(list)
    all_ops: dict[tuple[str, str], list[Overload]] = defaultdict(list)
    all_ns: set[str] = set()
    all_consts: dict[tuple[str, str], Constant] = {}
    all_enums: dict[tuple[str, str], Enum] = {}

    for rel in SOURCE_HEADERS:
        path = PROJECT_ROOT / rel
        print(f"Parsing {rel} ...", flush=True)
        tu, errs = parse_header(path, args)
        if errs:
            print(f"  WARNING: {len(errs)} error diagnostics:", file=sys.stderr)
            for e in errs[:5]:
                print(f"    {e.location}: {e.spelling}", file=sys.stderr)
            return 1
        types, funcs, ops, ns, consts, enums = collect(tu)

        for t in types:
            all_types.setdefault((t.namespace, t.name), t)
        for k, v in funcs.items():
            all_funcs[k].extend(v)
        for k, v in ops.items():
            all_ops[k].extend(v)
        all_ns |= ns
        for c in consts:
            all_consts.setdefault((c.namespace, c.name), c)
        for e in enums:
            all_enums.setdefault((e.namespace, e.name), e)

    manifest = Manifest(
        schema_version=SCHEMA_VERSION,
        generated_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        source_files=SOURCE_HEADERS,
        value_t="double",
        namespaces=sorted(n for n in all_ns if n.startswith("hd")),
        types=sorted(all_types.values(), key=lambda t: (t.namespace, t.name)),
        functions=to_groups(all_funcs),
        operators=to_groups(all_ops),
        constants=sorted(all_consts.values(), key=lambda c: (c.namespace, c.name)),
        enums=sorted(all_enums.values(), key=lambda e: (e.namespace, e.name)),
    )

    out_path = Path(args_ns.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(manifest.to_json(), encoding="utf-8")

    print(f"\nManifest written to: {out_path}")
    print(f"  types:      {len(manifest.types)}")
    print(
        f"  functions:  {len(manifest.functions)} unique names "
        f"({sum(len(g.overloads) for g in manifest.functions)} overloads)"
    )
    print(
        f"  operators:  {len(manifest.operators)} unique names "
        f"({sum(len(g.overloads) for g in manifest.operators)} overloads)"
    )
    print(f"  constants:  {len(manifest.constants)}")
    print(f"  enums:      {len(manifest.enums)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

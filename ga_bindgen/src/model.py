"""Manifest schema (dataclasses → JSON).

A manifest is the canonical, language-neutral description of the C++ API
surface that ga_bindgen has scanned from `ga/`. Emitters consume manifests;
the scanner produces them.

Schema version 1.
"""
from __future__ import annotations

import dataclasses as dc
import json
from typing import Any


SCHEMA_VERSION = 1


@dc.dataclass
class Field:
    name: str
    type: str  # textual type as written, e.g. "T" or "value_t"


@dc.dataclass
class Constructor:
    param_types: list[str]  # textual, e.g. ["T", "T", "T"]


@dc.dataclass
class TypeAlias:
    """A user-facing typedef like `vec3d -> Vec3_t<double, vec3d_tag>`."""
    namespace: str            # e.g. "hd::ga"
    name: str                 # e.g. "vec3d"
    underlying: str           # e.g. "Vec3d<value_t>"
    canonical_underlying: str # e.g. "hd::ga::Vec3_t<double, hd::ga::vec3d_tag>"
    fields: list[Field] = dc.field(default_factory=list)
    constructors: list[Constructor] = dc.field(default_factory=list)
    base_class: str | None = None  # short template name of base, e.g. "Vec2d"
                                   # for classes like `Point2d : public Vec2d<T>`
    source_file: str = ""
    source_line: int = 0


@dc.dataclass
class Overload:
    return_type: str
    param_types: list[str]
    source_file: str
    source_line: int


@dc.dataclass
class FunctionGroup:
    """All overloads sharing a fully-qualified name (e.g. `hd::ga::ega::gpr`)."""
    namespace: str
    name: str  # short name, no namespace
    overloads: list[Overload] = dc.field(default_factory=list)


@dc.dataclass
class Constant:
    """A const variable like `e1_2d` or `I_3d` declared at namespace scope.
    Source: ga_usr_consts.hpp."""
    namespace: str         # e.g. "hd::ga::ega"
    name: str              # e.g. "e1_2d"
    type: str              # e.g. "const Vec2d<value_t>"
    source_file: str = ""
    source_line: int = 0


@dc.dataclass
class Manifest:
    schema_version: int
    generated_at: str          # ISO timestamp
    source_files: list[str]    # relative to project root
    value_t: str               # "double" in v1
    namespaces: list[str]
    types: list[TypeAlias]
    functions: list[FunctionGroup]
    operators: list[FunctionGroup]
    constants: list[Constant] = dc.field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(dc.asdict(self), indent=2)

    @classmethod
    def from_json(cls, text: str) -> "Manifest":
        d: dict[str, Any] = json.loads(text)
        return cls(
            schema_version=d["schema_version"],
            generated_at=d["generated_at"],
            source_files=d["source_files"],
            value_t=d["value_t"],
            namespaces=d["namespaces"],
            types=[TypeAlias(
                namespace=t["namespace"], name=t["name"],
                underlying=t["underlying"],
                canonical_underlying=t["canonical_underlying"],
                fields=[Field(**f) for f in t["fields"]],
                constructors=[Constructor(**c) for c in t["constructors"]],
                base_class=t.get("base_class"),
                source_file=t.get("source_file", ""),
                source_line=t.get("source_line", 0),
            ) for t in d["types"]],
            functions=[_fg(g) for g in d["functions"]],
            operators=[_fg(g) for g in d["operators"]],
            constants=[Constant(**c) for c in d.get("constants", [])],
        )


def _fg(g: dict[str, Any]) -> FunctionGroup:
    return FunctionGroup(
        namespace=g["namespace"], name=g["name"],
        overloads=[Overload(**o) for o in g["overloads"]],
    )

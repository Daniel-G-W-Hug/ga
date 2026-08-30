-- test_grades.lua - gr() and rgr() across every algebra
--
-- gr(x) is the grade; rgr(x) = n - gr(x) is the anti-grade, the number of
-- dimensions the object does NOT populate. The interesting part is where each is
-- DEFINED, because that is a deliberate rule rather than an accident:
--
--   gr   -- every algebra, every graded type (32 overloads in C++)
--   rgr  -- pga2dp/3dp and cga2dc/3dc only (20 overloads)
--
-- The criterion for rgr is not the metric (anti-grade is just n - grade and needs
-- none) but the ANTI-PRODUCT family: those four algebras carry rgpr, rdot and
-- rcmt, and the anti-grade is their grade bookkeeping. ega2d/3d and sta4ds carry
-- only rwdg and define no rgr -- sta4ds did until 2026-08-30, with no caller and
-- no test in any language. The absence is asserted below so it cannot creep back.

print("=======================================================")
print("Testing gr() / rgr() across all algebras")
print("=======================================================")

-- {table, type, ctor arity, expected grade, algebra dimension or nil if no rgr}
local T = {
    {"ega", "scalar2d", 1, 0, nil}, {"ega", "vec2d", 2, 1, nil},
    {"ega", "pscalar2d", 1, 2, nil},
    {"ega", "scalar3d", 1, 0, nil}, {"ega", "vec3d", 3, 1, nil},
    {"ega", "bivec3d", 3, 2, nil}, {"ega", "pscalar3d", 1, 3, nil},

    {"pga", "scalar2dp", 1, 0, 3}, {"pga", "vec2dp", 3, 1, 3},
    {"pga", "bivec2dp", 3, 2, 3}, {"pga", "pscalar2dp", 1, 3, 3},
    {"pga", "scalar3dp", 1, 0, 4}, {"pga", "vec3dp", 4, 1, 4},
    {"pga", "bivec3dp", 6, 2, 4}, {"pga", "trivec3dp", 4, 3, 4},
    {"pga", "pscalar3dp", 1, 4, 4},

    {"sta", "scalar4ds", 1, 0, nil}, {"sta", "vec4ds", 4, 1, nil},
    {"sta", "bivec4ds", 6, 2, nil}, {"sta", "trivec4ds", 4, 3, nil},
    {"sta", "pscalar4ds", 1, 4, nil},

    {"cga", "scalar2dc", 1, 0, 4}, {"cga", "vec2dc", 4, 1, 4},
    {"cga", "bivec2dc", 6, 2, 4}, {"cga", "trivec2dc", 4, 3, 4},
    {"cga", "pscalar2dc", 1, 4, 4},
    {"cga", "scalar3dc", 1, 0, 5}, {"cga", "vec3dc", 5, 1, 5},
    {"cga", "bivec3dc", 10, 2, 5}, {"cga", "trivec3dc", 10, 3, 5},
    {"cga", "quadvec3dc", 5, 4, 5}, {"cga", "pscalar3dc", 1, 5, 5},
}

local function mk(name, arity)
    local ty = _G[name]
    assert(ty ~= nil, "type not bound in lua: " .. name)
    local a = {}
    for i = 1, arity do a[i] = i end
    return ty.new(table.unpack(a))
end

print("\n1. gr() returns the grade")
print("-------------------------")
local n_gr = 0
for _, e in ipairs(T) do
    local tbl, name, arity, want = e[1], e[2], e[3], e[4]
    assert(_G[tbl].gr ~= nil, tbl .. ".gr is not bound")
    local g = _G[tbl].gr(mk(name, arity))
    assert(g == want, string.format("%s.gr(%s) = %d, want %d", tbl, name, g, want))
    n_gr = n_gr + 1
end
print(string.format("  %d graded types, all correct", n_gr))

print("\n2. rgr() where defined: gr + rgr == dimension of the space")
print("----------------------------------------------------------")
local n_rgr = 0
for _, e in ipairs(T) do
    local tbl, name, arity, want, dim = e[1], e[2], e[3], e[4], e[5]
    if dim ~= nil then
        assert(_G[tbl].rgr ~= nil, tbl .. ".rgr is not bound")
        local v = mk(name, arity)
        local g, r = _G[tbl].gr(v), _G[tbl].rgr(v)
        assert(g + r == dim, string.format("%s.%s: gr=%d rgr=%d, sum %d, want %d",
                                           tbl, name, g, r, g + r, dim))
        n_rgr = n_rgr + 1
    end
end
print(string.format("  %d types in pga/cga, identity holds", n_rgr))

print("\n3. rgr() is ABSENT for ega and sta (only rwdg, no anti-product family)")
print("----------------------------------------------------------------------")
local n_abs = 0
for _, e in ipairs(T) do
    local tbl, name, arity, _, dim = e[1], e[2], e[3], e[4], e[5]
    if dim == nil then
        local v = mk(name, arity)
        local reachable = _G[tbl].rgr ~= nil
                          and pcall(function() return _G[tbl].rgr(v) end)
        assert(not reachable, tbl .. ".rgr(" .. name .. ") exists but should not")
        n_abs = n_abs + 1
    end
end
print(string.format("  %d ega/sta types, rgr correctly unavailable", n_abs))

print("\nOK - grades: " .. n_gr .. " gr, " .. n_rgr .. " rgr, " .. n_abs .. " absent")

-- test_operators.lua - the +/- type matrix, per algebra
-- Every ordered pair of an algebra's own types must support + and -, and the
-- result must be independent of how the operands are ordered. Both halves matter
-- and neither implies the other:
--
--   * REACHABILITY. Lua takes __add from the LEFT operand's metatable and falls
--     back to the right one only if the left has none, so an overload registered
--     in the wrong usertype block compiles, reads correctly and can never fire.
--     ga_lua reached 88 of 428 pairs before 2026-08-30 with nothing noticing.
--
--   * SELECTION. That a call succeeds says nothing about WHICH overload sol
--     picked. Addition commutes and subtraction anticommutes, and both are
--     type-stable here (the narrowest type holding the grade union is symmetric
--     in the operands), so a wrongly selected overload shows up as a value or
--     type mismatch in the identities below.

print("=======================================================")
print("Testing the +/- type matrix (all algebras)")
print("=======================================================")

local T = {}
T["ega2d"] = {"scalar2d", "vec2d", "pscalar2d", "mvec2d_e", "mvec2d"}
T["ega3d"] = {"scalar3d", "vec3d", "bivec3d", "pscalar3d", "mvec3d_e", "mvec3d_u", "mvec3d"}
T["pga2dp"] = {"scalar2dp", "vec2dp", "bivec2dp", "pscalar2dp", "dualnum2dp", "mvec2dp_e", "mvec2dp_u", "mvec2dp"}
T["pga3dp"] = {"scalar3dp", "vec3dp", "bivec3dp", "trivec3dp", "pscalar3dp", "dualnum3dp", "mvec3dp_e", "mvec3dp_u", "mvec3dp"}
T["sta4ds"] = {"scalar4ds", "vec4ds", "bivec4ds", "trivec4ds", "pscalar4ds", "mvec4ds_e", "mvec4ds_u", "mvec4ds"}
T["cga2dc"] = {"scalar2dc", "vec2dc", "bivec2dc", "trivec2dc", "pscalar2dc", "mvec2dc_e", "mvec2dc_u", "mvec2dc"}
T["cga3dc"] = {"scalar3dc", "vec3dc", "bivec3dc", "trivec3dc", "quadvec3dc", "pscalar3dc", "mvec3dc_e", "mvec3dc_u", "mvec3dc"}

local order = {"ega2d", "ega3d", "pga2dp", "pga3dp", "sta4ds", "cga2dc", "cga3dc"}

-- build a NON-ZERO instance: zeros would satisfy the identities trivially
local function mk(name, seed)
    local ty = _G[name]
    assert(ty ~= nil, "type not bound in lua: " .. name)
    for n = 1, 32 do
        local a = {}
        for i = 1, n do a[i] = seed + i * 0.25 end
        local ok, v = pcall(function() return ty.new(table.unpack(a)) end)
        if ok then return v end
    end
    error("no component constructor found for " .. name)
end

print("\n1. Reachability: every ordered pair supports + and -")
print("---------------------------------------------------")
local pairs_ok = 0
for _, alg in ipairs(order) do
    local ts = T[alg]
    for _, a in ipairs(ts) do
        for _, b in ipairs(ts) do
            local x, y = mk(a, 1.0), mk(b, 7.0)
            local ok_add = pcall(function() return x + y end)
            local ok_sub = pcall(function() return x - y end)
            assert(ok_add, alg .. ": " .. a .. " + " .. b .. " is not reachable")
            assert(ok_sub, alg .. ": " .. a .. " - " .. b .. " is not reachable")
            pairs_ok = pairs_ok + 1
        end
    end
    print(string.format("  %-8s %3d pairs", alg, #ts * #ts))
end
assert(pairs_ok == 428, "expected 428 ordered pairs, saw " .. pairs_ok)
print(string.format("  total    %3d pairs, all reachable", pairs_ok))

print("\n2. Selection: a + b == b + a and a - b == -(b - a)")
print("--------------------------------------------------")
local ident = 0
for _, alg in ipairs(order) do
    for _, a in ipairs(T[alg]) do
        for _, b in ipairs(T[alg]) do
            local x, y = mk(a, 1.0), mk(b, 7.0)
            assert(tostring(x + y) == tostring(y + x),
                   alg .. ": " .. a .. " + " .. b .. " does not commute")
            assert(tostring(x - y) == tostring(-(y - x)),
                   alg .. ": " .. a .. " - " .. b .. " does not anticommute")
            ident = ident + 2
        end
    end
end
print(string.format("  %d identities hold", ident))

-- Negative zero would break the string comparisons above; that it does not is
-- itself a check that every to_string goes through the C++ formatter, which
-- suppresses it (four of them did not until 2026-08-30).
print("\nOK - operator matrix: " .. pairs_ok .. " pairs, " .. ident .. " identities")

Please analyse the subfolder "ga_prdexp/src_prdexpr" and within that subfolder the four
files ga_prdxpr_ega2d.hpp, ga_prdxpr_ega3d.hpp, ga_prdxpr_pga2dp.hpp and
ga_prdxpr_pga3dp.hpp. For our purpose today only the prd_rules*variables are relevant. Don't condsider the mvec_coeff that are created after the prd_rules* for now. The only mvec_coeff that is relevant is the first one in each file, because is determines the basis elements of the algebras to be generated. The relevant maps, that are provided by the user up-to-now, determine the geometric algebras that are created in the GA library. Since I am starting to create a new algebra, I want to avoid the tedious process and generate those maps automatically, based on some calculation rules. The maps contain the products of all basis elements with each other and a simplification step. The key rules for the vector basis elements are given in the comments at the beginning of the file and in addition there is the rule that basis elements multiplied by a scalar (prepresented by "1") are the basis elements themselves. The sign is determined by how many index swaps are need to get back to a permissible basis element. Odd number of swaps yield a minus, even number of swaps yield a plus. Do you think you could provide routines where, given the inputs "basis coefficients" and the "rules for multiplication of the basis vectors" (now in comments), the maps that are now user provided,
can be created by calling these routines? What questions are still open and what would be
the proposed steps to try that for one algebra like ega2d as a first step? The other
algebras should be created later on. The given user inputs in the four files mentioned
above shall of course be used as a reference. Please be aware that the basis
multiplication rules so far contain only results of "+1" or "0", but currently I am
preparing space-time algebra which will contain also "-1" as possible result (see
ga_prdxpr_sta3d.hpp, but only worked on until line 105).

I would like you to make a refinement: There are some definitions in ga_prdxpr_common.hpp (and maybe also useful things in ga_prdxpr_common.cpp, but I am not sure). Please use the ga_prdxpr_common.hpp header and the definitions contained in there in order to get consistent usage of string definitions. Thanks!

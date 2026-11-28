# Hint for application

- **Single include file for the user** (depending on the GA needed):
  - ga/ga_ega.hpp **or**
  - ga/ga_pga.hpp

- The value type `value_t` is **double by default; float is selected with the compile
  definition `_HD_GA_VALUE_T_FLOAT`**, e.g. `target_compile_definitions(<target> PRIVATE
  _HD_GA_VALUE_T_FLOAT)`. Set it for every translation unit of the program, or for none
  (see ga/ga_value_t.hpp); the build checks the float configuration with its own test
  targets

- Just for easy user reference following includes are provided at the top level:
  - ga/ga_usr_consts.hpp
  - ga/ga_usr_types.hpp
  - ga/ga_usr_utilities.hpp

    **These three include files are already included automatically by ga/ga_ega.hpp or
    ga/ga_pga.hpp respectively, and thus must not be included by the user again.**

    All other include files are not intended for direct use by the user of the ga library.

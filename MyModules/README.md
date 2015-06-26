# MyModules

"MyModules" directory is aimed at managing all process modules which will be included while compilation of a THMC simulator. CMake will search modules in sub-directories and all modulues are included in a compilation by default. To exclude specific modules from the compilation, please modify InactiveModules.cmake in this directory.

Each module directory should cotain the followings
- ConfigureModule.cmake: register process, set up source and test codes
- source dir: include source codes
- tests dir: include test codes based on google test framework

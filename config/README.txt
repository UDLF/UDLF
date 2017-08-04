This directory stores the internal configuration files for validation!
These files specify the parameters that the user can change in the config.ini file.
WARNING: Do not change anything here unless you really know what you are doing!

Notice that:
1) The *.conf files are internal configuration files that are embedded into the binary at compile time.
2) All the files must follow the structure R"================(conf file content goes here)================"
3) While the "general.conf" specifies the global parameters (required by all the methods), the other files specify method specific parameters.
4) Method specific parameters must start with "PARAM_METHODNAME_"
5) The files are loaded into char buffers by Core/Conf.cpp and are used in the Core/Validation class.

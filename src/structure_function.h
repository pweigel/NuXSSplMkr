#ifndef __STRUCTURE_FUNCTION_H
#define __STRUCTURE_FUNCTION_H

namespace nuxssplmkr {
    struct SFConfig {
        double test;
    };

    class StructureFunction {
        private:

        public:
            StructureFunction(SFConfig sf_config);
            ~StructureFunction();
    };
}

#endif
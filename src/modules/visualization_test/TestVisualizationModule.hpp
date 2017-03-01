/*
 * Visualization module
 */

#ifndef ALLPIX_TEST_VISUALIZATION_MODULE_H
#define ALLPIX_TEST_VISUALIZATION_MODULE_H

#include <memory>
#include <string>

#include "../../core/module/Module.hpp"
#include "../../core/config/Configuration.hpp"

class G4UIsession;
class G4VisManager;

namespace allpix{
    // define the module to inherit from the module base class
    class TestVisualizationModule : public Module{
    public:
        // provide a static const variable of type string (required!)
        static const std::string name;
        
        // constructor should take a pointer to AllPix, a ModuleIdentifier and a Configuration as input
        TestVisualizationModule(AllPix *apx, ModuleIdentifier id, Configuration config);
        ~TestVisualizationModule();
        
        // method that will be run where the module should do its computations and possibly dispatch their results as a message
        void init();
        void run();
        
    private:
        // configuration for this module
        Configuration config_;
        
        // Shared pointer to the session
        std::shared_ptr<G4UIsession> session_g4_;
        
        std::shared_ptr<G4VisManager> vis_manager_g4_;
    };
}

#endif /* ALLPIX_TEST_VISUALIZATION_MODULE_H */

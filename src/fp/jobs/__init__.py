import importlib
import pkgutil
import os
import sys

# __all__ will contain all imported names from the submodules
__all__ = []

# Get the path of the current package
package_path = os.path.dirname(__file__)

# Iterate through all modules in the package
for module_info in pkgutil.walk_packages(path=[package_path], prefix=__name__ + "."):
    # Import the module
    module = importlib.import_module(module_info.name)
    
    # Get the module's public variables (i.e., those without a leading underscore)
    module_vars = {k: v for k, v in vars(module).items() if not k.startswith('_')}
    
    # Update the globals of the current __init__.py module with the variables from the imported module
    globals().update(module_vars)
    
    # Add these names to __all__ for explicit export
    __all__.extend(module_vars.keys())

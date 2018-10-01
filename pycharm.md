To run the mantid_total_scattering script with PyCharm, you will need to:
1. Setup the project configuration
1. Setup the test/run configuration

### Project configuration
* Go to project settings (`File > Settings` or `Ctrl+Alt+S`)
* In Project Interpreter add Mantid Python:
  * Select show all from the drop down menu
  * Click the `+` button
  * Click the `...` for the `Base interepreter` and navigate to `git/mantid/external/src/ThirdParty/lib/python2.7/python.exe` 
* In Project Structure add the mantid source:
  * Click `+ Add Content Root`
  * Navigate to `build/mantid/bin/Debug`
  * Click `Apply`


### Test configuration
* Create a new configuration (using the configurations drop down menu in the top right and clicking `Edit Configurations...`) 
* Click `+` in the top left
* `Python tests > pytest`
  * Name the test something sensible in the top box
  * Check `Script path`
  * Click the folder symbol in the right on the input field and navigate to the test you want to run
  * Add the follwoing environment variables (ensure to keep the semi colons)
    * PATH : 
      * `<path-to-build>/bin/Debug;`
      * `<path-to-git>/external/src/ThirdParty/bin;`
      * `<path-to-git>/external/src/ThirdParty/lib/qt4/bin;`
      * `<path-to-git>/external/src/ThirdParty/lib/qt5/bin;%PATH%`
  * Ensure that the `Python interpreter` is the mantid interpreter we define in `Project Configuration` above.

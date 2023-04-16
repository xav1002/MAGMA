%% Development Notes

% Component Name Notes: Subject_Action_Component (e.g. SpecAB = Species Add
%   Button)
% Spec = Species
% Chem = Chemical
% Environ = Environment
% Cust = Custom
% F = text Field
% DD = Drop Down
% B = Button
% A = Add
% U = Update
% R = Remove
% T = Table
% Tab = Tab
% Tree = Tree
% C = Checkbox
% L = Label
% P = Panel
% Subst = Substrate

% SpecT:
% Col1: spec number
% Col2: spec name

% Model Function Variable Nomenclature
% EX: K12
% Interpretation: Half-Saturation of first Chemical for 

% ODESys.updateSpecies:
%   field: string that corresponds to a property of Species to update

% 11/11/22:
% The ODESys class needs to be able to store all of the Species and
% Chemicals. Update this class to have a property storing the Species and
% Chemicals, and methods to edit that property. The editing methods should
% then chain to update the dydt (system functions property) of the ODESys
% class, so there's always an up-to-date ODE system.
% All unit conversions for inputs will be done in the callbacks for each
% UI component. All unit conversions for outputs will be done in the
% graphing stage. Edit Species: callback => func within ODESys => update
% within Species => update within func within ODESys

% 11/13/22:
% Need to figure out how to import classes in MATLAB: classes are not on
% the correct path (possibly because files are stored in OneDrive? Try
% moving to local directory).
% Write code to add/store/delete species, maintain correct order in table
% and in the ODESys storage.

% 11/14/22:
% dictionary data type doesn't work, use containers.Map
% Successfully added species to SpecT (SpecABPushed())

% 11/14/22:
% containers.Map and dictionaries can't hold custom classes
% Use struct data type to hold Species and Chemical classes
% Need to replace spaces in species name with underscores to use as struct
% data type field name

% 11/16/22:
% Working on removing species button
% Made ODESys into a handle object (default is a value object, see
% documentation) and this enabled persistence of changes to ODESys
% properties
%   Also made Species and Chemical classes into handle classes
% Created inputs for Environment, will allow ODESys to store default and
% custom environments
%   Created Environment class

% 11/17/22:
% Creating default values for parameters
%   Created SpecDefaults class that is instantiated in ODESys constructor
%   and has properties that store all of the default parameters
% Implementing add/update/remove functions in ODESys for Environment class
% Created EnvDefaults class that is instantiated in ODESys constructor and
% has properties that store all default params for Environment
% Created "Reset to Default" and "Set as Default" buttons for Environments

% 11/20/22:
% Coding unit conversions for Environmental parameter updates.
% Finished implementing removeSpecies functionality button.
% Implementing Environment update functionality.
% Active Environment name is stored in ODESys.activeEnv property
%   Use this property to index the ODESys.environs struct when pulling data
%   for plotting.
% Finished code in ODESys for updateEnvironment.
% Started code to update the display tables for Environments.

% 11/21/22:
% Finished code to display tables for Custom Environment default values.
% Starting code to display tables for active Environment values.
% Updated ODESys.getDefaultEnvironParams() to return table-ready values,
% don't need to format the data on the fly.
% Finished Environment tables updating functionalities.

% 11/23/22:
% When Environment class is instantiated, the lightFunc and tempFunc are
% created as function_handle; the data for the Custom Environment is not
% yet instantiated into a class, so it comes out just as doubles
% Fixed issue for displaying function_handle in table (func2str)/
% Finished first run of development for Species and Environment tab.
%   Started majority of groundwork for rest of application.
% Started looking into popups (multi-window app) for chemical model
% specifications. Created new sub-app, Chem_Model_Popup.mlapp.

% 11/26/22:
% In Chemical model specification window:
%   1. Field to enter custom model equation;
%   2. Dropdown to select model;
%   3. Field to display model equation in fancy text;
%   4. Field(s) to enter parameter names;
%   5. Table to display 
% Will display in table in Chemical tab: name of model chosen for each
% species for each interaction functionality
% Updating the Chemical values in ChemT will have data go through Chemical
% class in ODESys before being displayed in table. This is because it is
% more complex to display the model data for chemicals.
% Brainstormed how to create Chemicals UI:
%   Need to be able to display each Chemical, all the Species it interacts
%   with, how the Chemical interacts with each Species, and the model used
%   to represent that interaction.
%   Split into 2 tabs? Add Chemicals + Specify Models
%   First tab can just be adding Chemicals + init concentration (table
%   shows all chemicals)
%   Second tab can be go through chemicals one-by-one, tables show data one
%   Chemical at a time.

% 11/27/22:
% Created new tab for Chemical interaction model specification.
% Started implementing adding/removing Chemical functionality.
% Using replace() for replacing '^' with '_', regexprep sometimes doesn't
% work.

% 11/30/22:
% Implementing Interaction Model function editing.
%   ModelAB will run ODESys.updateModel(). This function takes input of
%   Model Name, Chemical, Interaction Type, Species, and Function. The
%   Parameter names will be pulled from a temporary store (see below).
%   tempParamStore will be the cell array that stores the individual
%   parameter name string. These are added to tempParamStore when the
%   ModelParamAB is pushed, and removed when the Chemical, Interaction, or
%   Species is changed. This store is persisted when the ModelAB is pushed.
% Implementing model setting and storage functions in Species, Chemical,
% and ODESys classes.
%   In ODESys, still need to put in the storage for parameter name strings
%   for preset model functions.
%   Implemented str2func in ODESys with parameter name strings.

% 12/1/22:
% Adding functionality to return (from ODESys) all of valid param names for
% the Interaction Models with getModelParamNames(). getModelParamNames()
% will be run from the application side, and directly pass data into app.ModelParamT
% In the process of adding parameters from Environment to the Interaction
% Model Parameter Table

% 12/2/22:
% Writing code to set model function in Species and Chemical from UI side.
%   This is called by the ModelAB being pushed.
% Writing code to set model function in the ModelF.
%   This is called by changes to the ModelDD value.
%   When this is not on "Custom Model", it will pull variables names based
%   on the Chemical and Species that are being used.
% Writing code for LaTeX display of model functions.
% Should use normal function for ODESys composite function?
%   Or should use anonymous function?
%   Still need to separate variables for LaTeX display.
%   Don't need to create full anonymous function handle for each individual
%   model function.
%   Working on ODESys.recursiveFrac to recursively create fraction display
%   for LaTeX function display.

% 12/10/22:
% Using latex() function to turn function input into LaTeX output.
% Updating Linear model.
% Updating the "Model Function" panel to be more clear about how each
% individual function is used in the composite growth function for each
% Species and Chemical.
% Creating function in ODESys to return functions to UI (getModelFunc).
% Default values in Interactions Models Panel
% 2 portions of functions: "multiplied to mu,max" and a "added to mu,max".
%   Currently, only the "multiplied to mu,max" is able to be added.
% Created a G function input (add to) portion, will simultaneously update
% with F portion (multiply to) portion.

% 12/24/22:
% Creating "Home" tab to give user option to choose from existing model
% scheme or to create new custom scheme.
%   User can still edit existing scheme if they choose to use it.
% Wrote recursive function for enabling or disabling all child components
% of each tab.
% Created model summary tab for user to view final model.
%   Depending on if user chooses to use existing model or create custom
%   model, the app will take the user to the correct tab.
% Can create pop-ups for instructions with the Panel component.

% 1/3/23:
% How can the users store their new custom model schemes? Is there a
% backend system that can be used for this? Maybe a github repo?
%   Can integrate external files into MaGKMA.mlapp with input arguments
%   into the startupfunction for the application.
% Design the app as if this will work, when we get to the point where we
% need to implement this, we can get someone from CS to help us.
% Created general set up for Model Scheme management on Home Page.
% Currently working on saving and assembling ODE models in Species and
% Chemical classes.
%   Including funcType ("F" or "G").
%   Conversion of function strings to function handles will happen when all
%   components of the entire function are collected.
% Also updating the model summary tables with ODESys.updateModelSummaryT().
% Upgraded ODESys.getModelFunc() to also return a G function value.
% Next Time: Write setModel() for Species and Chemical classes, then write
% code for assembling final function_handle's and display them on Model
% Summary tab.

% 1/4/23:
% When creating the aggregate functions for dydt, can only include the y
% (values of the dependent variables) in the parameters that are passed
% in to dydt function.
%   Need to have system to correspond each Species and Chemical
%   concentration to a given index in y(t).
%   Need to have system to set values for growth parameters in ODESys.
%       Can maintain these parameters in Species and Chemical classes ass
%       needed.
%       System will be to store new models, variable, and parameter
%       names/values in each Species and Chemical class, then rely on
%       Species and Chemical number to maintain the index of each dependent
%       variable in y(t).
%   May need to pass in parameters (constants) into the dydt function
%   itself (pass an array of constants as one of the input parameters?).
%       Can do this, see syntax on "Pass Extra Parameters to ODE Function"
%       section in documentation for ODE45.
%   Need to figure out how to convert param names into p(index).
%   Will be passing raw user input into Species and Chemical classes to
%   convert function to terms of t, y, and p to send back to ODESys class.

% 1/6/23:
% The parameter numbers should correspond to the parameter numbers on the
% Growth Parameters tab table.
%   Adding Growth Parameters to array of params should be done in
%   ODESys.updateModel().
% Only dependent variables in the system should be displayed as inputs to
% each F and G function.
% Updates from Species and Chemical classes to ODESys (of params or
% functions) should create a new array of all params or functions from
% scratch each time, not trying to keep track of and edit existing array.
% Updated ODESys.updateModel() to only output dependent variables (or time,
% IV) as inputs to F and G functions.
% Finished code to input dependent varNames, paramNames (string arrays)
% and function strings to Species.setModel() and Chemical.setModel().
% All of the parameters of the system need to be consolidated into a 1xn
% array. This way, the parameters can all be passed in to the ode45
% function as a single array. The parameters will be organized by Speices
% (e.g. for each Species, the list of parameters displayed in the Growth
% Parameters table will haev those corresponding numbers in the appropriate
% Species class parameter array). The parameters will need to have their
% corresponding numbers be coordinated in the Chemical classes.
% Design for keeping track of parameters across Species and Chemicals:
%   Add the function and params to Species.
%   Species.setModel() should return the numbers for the parameters in
%   the specified Species.
%   The Chemical class can then track the parameter number and the Species
%   number (which it receives as a parameter to Chemical.setModel()) to
%   calculate the correct number for the parameter when the ode45 function
%   is assembled.

% 1/7/23:
% Work on storing the functions in each Species and Chemical class.
% When setting function in Species class, will create a parameter struct
% for every param name.
% When the user sets a relationship between a Species and Substrate to
% "Non-Interactive With", that essentially means the function that was
% existing between the Species and Substrate needs to be deleted.
%   This means that parameter also needs to be deleted.
% For function default parameter values, units, names, etc, store/access
% these data from the FuncDefaults.m file.
% Working on updating the Growth Parameters table to have the parameters
% corresponding to the selected Species.
% Finished first run of updating the parameters in the growth parameter
% table. Need to debug.

% 1/8/23:
% Debug Growth Parameter table.
% Look into MongoDB.
%   MongoDB Atlas is a cloud based DB service that can be run through AWS,
%   Azure, or Google Cloud. Can be free for a small amount of storage, but
%   can also be paid for for a larger amount of storage and support. I can
%   probably learn how to use this, and Garrett probably has a better
%   handle than I do.
%   Look into MongoDB Atlas documentation on how to set up database and for
%   pricing options.
%   There is Matlab documentation on using MongoDB with MATLAB in C++
%   interface.
% The Dropdown component can be used as a type of dynamic search bar.
%   Will need to be careful with debugging.

% 1/17/23:
% Looked into webread() function:
%   Can pull data from the Firebase DB, and can custom query in the
%   webread() function (may be faster) or can pull entire database and do
%   filtering on the MATLAB side.
% Need from Garrett:
%   1. What URL should I send POST requests to?
%   2. What URL should I send GET requests to?
%   3. Is there a certain format that I should be using?
%   4. Are there any DON'Ts regarding database security or data usage?
% Method for storing the parameters and have them correspond to the
% functions:
%   FuncParam Cell Array:
%   Each cell contains a struct with fields spec, chem, func and params, will be
%   property of ODESys.
%   func is a string of the function (without the parameters being replaced
%   with p(i) yet
%   params is a cell array of parameter structs, containing param num,
%   name, sym, val, units
%   spec and chem specifies the spec and chem
%   Need to create functions to update the numbers for parameters every
%   time there is a change to the functions.
%   The info path is MAGMA.mlapp => ODESys => Species (or Chemical) =>
%   ODESys => Species (or Chemical)
%       Initial user input => Pass to ODESys (updateModel()) =>
%       Species.setModel() => return value to ODESys to update FuncParam
%       cell array => pass values into spec/chem when generating functions.
% const firebaseConfig = {
%   apiKey: "AIzaSyC9BPzYSHZADIl-ihS8mrk6Hcw7gIHPOmk",
%   authDomain: "impact-db.firebaseapp.com",
%   projectId: "impact-db",
%   storageBucket: "impact-db.appspot.com",
%   messagingSenderId: "625778730034",
%   appId: "1:625778730034:web:ded85e74593cc3e554482b",
%   measurementId: "G-MQJE957JGM"
% };
% Look up tutorials
% Cloud functions

% 1/20/23:
% There seems to be issue with using webread and webwrite to access
% FireStore. An alternative is to use a NodeJS script to do the
% communication with Firestore, then use Matlab to run those NodeJS
% scripts.

% 1/23/23:
% Shouldn't use a FuncParam cell array in ODESys, should store functions in
% Spec and Chem classes. This will allow easier handling of duplicate
% parameters. In order to handle parameters that appear for chem only and
% not spec, should add Chem objects to Growth Parameter selection DD. Each
% parameter should have a local number (local to Spec or Chem) and global
% number (across all Spec/Chem). Creating a funcParams cell array for each
% Spec/Chem class.

% 1/24/23:
% Working on funcParam object handling functions in Species. Running into
% issues with get function. Maybe just get by each individual field, then
% use Matlab functions to cross compare each result?

% 1/25/23:
% Working on filtering/finding in cell arrays for handling params and funcs
% in funcObjs in Species.
% In function input, may want break down the functions into the
% relationship between 1 spec and 1 chem at a time, then use MATLAB to
% create more succinct overall functions once all of the independent
% relationships have been established.

% 1/26/23:
% Rewriting getGrthParam() in Species to return growth parameter data to
% ODESys for display. Need to properly write the funcParam creation code,
% to prevent empty structs from being thrown into the funcParams property.

% 1/27/23:
% Can use symbolic functions (use sym(function_handle) to convert anonymous
% func string to symbolic).
%   Can replace substitute symbols with different symbols with subs()
%   If the parameters can be made into variables, putting them into a
%   function_handle and then converting to symbolic will automatically put
%   the parameter values in
%   Can turn char into symbolic variable with syms 'char'
%   Then can use subs() to replace symbol with numerical value
%   To use a symbolic function with ode45, need to convert back to a
%   function_handle
%   Alternatively, keep everything in function_handle, and just replace the
%   text
%   Symbolic equations can be analytically differentiated with diff()
%   Use matlabFunction(symFunc) to convert symbolic eq to function_handle
%   Can convert string input into a symbolic function via str2sym()
% Workflow for user string input to final numeric function
%   User string => str2sym() => symbolic function with variables + params
%   all symbols => subs() => symbolic function with only variables left as
%   symbols and params subbed out for numerical values => simplify() =>
%   simplified symbolic function => matlabFunction() => numerical function
%   (function_handle)
%   Use symar() to find all symbolic variables in symbolic function, then
%   match them to the system variables, every sym that isn't a system
%   variable is therefore a parameter
% When pulling funcParam data to create functions, need to convert units of
% parameters to standard units
% Working on updating parameter value and unit from ODESys.
% Working on editing paramDD values.

% 1/31/23:
% Working on creating a demo video and getting the functionality to fully
% generate a modeling scheme.
% Successfully got the updating functionality working for Growth Parameters
% Table. There are bugs, but these can be fixed, namely that updating the
% model function causes the table to not clear previous params and only add
% to the existing params.

% 2/2/23:
% The reason the updateModel() function is being called twice in Species is
% because the call is made twice in ODESys, one for F func and one for G.
% Need to figure out how I want to structure the input of different types
% of functions.

% 2/3/23:
% Adjusting the tables that record the interactions with each species for
% each chemical will require storage of data in the Chemical classes.
% Need to decide how to do the function input:
%   1. More structured: FFunc and GFunc
%   2. Less structured: each Species-Chemical pair has a fully defined
%   interaction
%       Will require special logic to simplify/get the shared
%       (multiplication) parts of the functions correct.
%   3. Keep structured scheme, but only make the functions multiplied by
%   Xi and Ci: allows for the F and G structure, but gives a little more
%   flexibility.
%   4. Have different segments of input for different parts of the function
%       Maybe still F and G, but also have a dropdown to specify what type
%       of func is being inputted?
%   5. Each biological and chemical species have the functions defined on
%   its own, but allow MAGMA to track what other functions need to be
%   created based on previously created functions.
%       Create a status bar like HYSYS?
%       Have general schemes? e.g. direct linear relationship between
%       biological and chemical species change rates.
%           This is like a "Fluids package" in HYSYS.
%       Customizable for each relationship between biological species and
%       chemical species?
%       Need that flexibility, because the models are not always going to
%       be in the style that Shawn has.
%           Need to give the user options on the scheme style.

% Tang: Make LaTeX display dynamic on Species Interaction Models page.
% Todo for demo:
%   1. Fix function inputs
%       For now, go with F and G funcs with scaling done only with Xi, not
%       mu,max
%   2. Make LaTeX display dynamic
%   3. Display of final function on Model Summary Tab
%   4. Fix parameter table not properly clearing
%
% Demo Script:
% Hello, welcome to a demo for the Micro-Algae Growth Modeling Application,
% nicknamed MAGMA. The purpose of this application is to create a platform
% for the efficient creation of kinetic growth models for micro-algae
% cultivation systems that is easy to learn and use without the
% requirement of extensive MATLAB knowledge. The goal is to allow users to
% either create their own modeling schemes, or to make modifications to an
% existing scheme created by other users. Creating a model from scratch
% allows for vast flexibility and specificity in the model, while modifying
% an existing model gives quick access to a defined model structure and
% parameters. This demo will give a quick overview of how to build a model
% from scratch.
%
% Along the top of the application window, there will appear a set of tabs.
% The first tab is the Home Tab, where the functionality to load existing
% modeling schemes will eventually be implemented. The Environment Tab is the
% beginning of the model building process. It allows the user to define the
% environmental conditions for the model, such as incident light and
% temperature. Modifying these conditions is done by selecting a condition
% from the list, inputting the desired value and units, and clicking
% update. When applicable, environmental conditions can be defined as a
% function of time.
%
% The Species and Chemicals Tab allows the user to specify which biological
% and chemical species are involved in the model, along with their initial
% concentrations. The chemical species include any nutrients, growth
% inhibitors, and metabolic products. The biological and chemical species
% are specified by selecting a species from the preset list and setting the
% initial concentration. If a species is not in the preset list, it can be
% inputted as a custom species.
%
% The Interaction Models Tab allows the user to specify the individual
% relationship a chemical has with a biological species. The user will define the 
% chemical substrate, the relationship the substrate has with a biological species, 
% and the biological species. The user is then able to build a function
% representing the individual rate of change relationship of the concentration for the
% biological and chemical species selected. The functions built by the user
% are in the form of an F and a G function. For a given species, all of the F and G functions will be compiled 
% according to these formulae to create the full rate of change equation,
% and same for a given chemical.
% These individual relationship functions will ultimately have units of 1/time, as they
% will be scaled to concentration per time as shown here.
% The purpose of this method is to provide the user with flexibility in
% modeling, while maintaining an organized structure.
% For example, a basic Monod relationship can be established between
% nitrate and 2973. Clicking the "Update" button will save this function as
% the function representing the relationship of nitrate being a nutrient
% for 2973. This function utilizes the nitrate
% concentration, a dependent variable of the system, that is listed in the
% System Variable table here.
%
% The Growth Parameters Tab allows the user to specify the values and units
% for any parameters that are used in the rate of change functions that
% aren't independent or dependent variables of the system. For example, the
% half-saturation constant for nitrate for 2973 can be set to a value of
% 0.15 g/L.
%
% When all relevant individual relationships and parameters have been specified, MAGMA will compile these
% relationships into an aggregate ODE system to represent the model. The
% ODEs in this system will be displayed on the Model Summary Tab. This
% portion of the model is still in development, but it will look something
% like this. The user will then be able to view the compiled functions, and
% go back to previous Tabs to edit the relationships if necessary.
%
% The user can then go to the Plot Customization Tab to specify which plots they want to
% create. This part of the user interface has not yet been created, but it
% will give the user options for plots such as cell concentration over
% time, nutrient concentration over time, final cell concentrations, etc.
% After selecting and customizing the plotting output, the user can then
% view the plots on the Plot Display Tab.
% Thank you.
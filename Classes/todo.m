%% Todo
%% Species and Environment
% 1. Create option to change order of species in Species and Environment
% 2. Limit the characters allowed in species names
% 3. Function input for light, temperature variation over time
%   Use str2func()
% 4. Finish writing the Environments class
% 5. Create Update, and Remove functions in ODESys for Environments
% 6. Write default conditions and parameters
% 7. Write setDefault() functions for Species and Chemicals
% 8. Set Environment as active based on dropdown
% 9. Implement all dropdowns and buttons on Species and Environments tab
% 10. Check unit conversions
% 11. Figure out table updating with Custom Environments
% 12. Reset Gen Param to Default button functionality
% 13. Fix data types for EnvDefaults structs
% 14. Fix aesthetics: data types have different alignments on tables
%   Could also just change the table alignment settings
% 15. Write default values for Species
% 16. Fix unit conversion logic to only convert after table display is
% finished.

%% Interaction Model Tab
% 1. Set up UI for chemical model specifications.
%   Need to be able to display each Chemical, all the Species it interacts
%   with, how the Chemical interacts with each Species, and the model used
%   to represent that interaction.
%   Split into 2 tabs? Add Chemicals + Specify Models
%   First tab can just be adding Chemicals + init concentration + specifying which Species each Chemical interacts with (table
%   shows all chemicals), lists interaction type
%   Second tab can be go through chemicals one-by-one, tables show data one
%   Chemical at a time, allows user to specify model equations.
%   uilabel component can output LaTeX with Label.Interpreter = 'latex'
%   DD for saved models, Latex display for models, Input field for
%   parameter variable names, table to display variables, 
% 3. Create add/remove functionality for Chemicals
% 5. Write setModel() functions in Species and Chemical classes
% 6. Create params store for preset model functions
% 8. Write code to save the functions in Species and Chemicals
% 9. Fix bugs with function display
% 10. Create dynamic functionality for the Spec and SubstMomdelEXL and
% ModelFuncP title.
% 11. Write functions to dynamically update the model record tables, switch
% LaTeX display for functions based what Chemical and Species are selected.
% 12. Write functions to create tables for Growth Parameter names.
% 13. 

%% Growth Parameters Tab
% 1. Create table columns:
%   Number, Parameter Name, Equation Symbol, Value, Unit
% 2. Optimize by converting chars to strings whenever possible?
%   Beware, all Items for dropdowns are chars.

%% New Features:
% 1. Optimization of yield, titer, production rate, etc.
% 2. Cost estimates?
% 3. Piecewise functions

%% Model Template Upload Interface
% Functions:
%   1. Search Model Schemes
%       Dynamic search bar;
%       Alphabetical List;
%   2. Select Model Scheme from list
%   3. Add Custom Model Schemes
%   4. Display basic information about selected Model Scheme
%   5. Set Model Scheme
%   6. View Model Functions (button to Model Summary Tab)
%   7. 


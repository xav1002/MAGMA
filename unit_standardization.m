function val = unit_standardization(val,unit)
    if strcmp(unit,"g/L")
        % concentration
    elseif strcmp(unit,"mg/L")
        val = val./1000;
    elseif strcmp(unit,"mg/g")
        val = val./1000;
    elseif strcmp(unit,"g/mg")
        val = val.*1000;
    elseif strcmp(unit,"cells/mg")
        val = val.*1000;
    elseif strcmp(unit,"mg/cells")
        val = val./1000;
    elseif strcmp(unit,"hr")
        % time
    elseif strcmp(unit,"min")
        val = val ./ 60;
    elseif strcmp(unit,"sec")
        val = val ./ (60.*60);
    elseif strcmp(unit,"day")
        val = val.* 24;
    elseif strcmp(unit,"K")
        % temperature
    elseif strcmp(unit,"C")
        val = val + 273.15;
    elseif strcmp(unit,"F")
        val = (val - 32) * (5/9) + 273.15;
    elseif strcmp(unit,"R")
        val = ((val) * (5/9));
    elseif strcmp(unit,'L')
        % volume
    elseif strcmp(unit,'m^3')
        val = val .* 1E3;
    elseif strcmp(unit,'mL')
        val = val .* 1E-3;
    elseif strcmp(unit,"bar")
        % pressure
    elseif strcmp(unit,"kPa")
        val = val ./ 100;
    elseif strcmp(unit,"MPa")
        val = val .* 10;
    elseif strcmp(unit,'L/s')
        % volume per time
    elseif strcmp(unit,'m^3/s')
        val = val .* 1000;
    elseif strcmp(unit,'mL/s')
        val = val ./ 1000;
    elseif strcmp(unit,'bar/s')
        % pressure per time
    elseif strcmp(unit,'kPa/s')
        val = val ./ 100;
    elseif strcmp(unit,'MPa/s')
        val = val .* 10;
    elseif strcmp(unit,"m")
        % length
    elseif strcmp(unit,"cm")
        val = val ./ 100;
    elseif strcmp(unit,"mm")
        val = val ./ 1000;
    elseif strcmp(unit,"rpm")
        % frequency
    elseif strcmp(unit,"rad/s")
        val = val ./ (2.*pi.*60);
    end
end
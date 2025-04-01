% HFCONDPROP Initialize a set of high-frequency dielectric properties of ice for use with HFCOND.
%   
%   HFCP = HFCONDPROP returns the structure HFCP containing a default set
%   of values and uncertainties for the dielectric properties of meteoric
%   polar ice at high frequency (HF). This structure is a necessary input
%   to HFCOND, which calculates the HF conductivity of ice based on HFCP,
%   temperature, impurity concentrations and density. Table 1 of MacGregor
%   et al. [2007] lists the default values used, which are directly
%   calculated by HFCONDPROP based on the data compilation in that table.
%   
%   HFCP = HFCONDPROP(MODEL) assigns HFCP to a given model:
% 
%       'M07-std'       Mean compilation of dielctric properties in Table 
%                       1 of MacGregor et al. [2007] (default).
% 
%       'M07-wt'        Same as 'M07-std' except that values are calculated by
%                       weighting using their uncertainties following
%                       Bevington [1969] This model is not recommended, as
%                       per the reasons given by MacGregor et al. (2007),
%                       section 2.2.2.
% 
%       'M07-adj'       Same as 'M07-std' except that pure ice values are
%                       reassigned to values determined by Johari and
%                       Charette [1975], following the discussion in
%                       MacGregor et al. [2007]. This model generally
%                       results in higher predicted HF conductivities.
% 
%       'M07-adj-NH4'   Same as 'M07-adj' except that values for NH4+ are
%                       set to match 'W97'.
% 
%       'W97'           Model based Wolff et al. [1997] empirical DEP-chemistry
%                       relationship for the GRIP ice core that
%                       incorporates physical insights into HF conduction
%                       from Stillman et al. [2013].
% 
%       'W97-adj'       Same as 'W97' except corrects all values so that
%                       radar-inferred and borehole temperatures match in
%                       Greenland.
%   
%   HFCP = HFCONDPROP(MODEL,'PARAM1',VALUE1,'PARAM2',VALUE2,...) reassigns
%   PARAM1 to VALUE1, PARAM2 to VALUE2, etc.
%   
%   This function follows Table 1 and section 2 of:
%   
%   MacGregor, J.A., D.P. Winebrenner, H. Conway, K. Matsuoka, P.A.
%   Mayewski and G.D. Clow, 2007, Modeling englacial radar attenuation at
%   Siple Dome, West Antarctica, using ice chemistry and temperature data,
%   J. Geophys. Res., 112, F02008
%   
%   Additional references:
%   
%   Johari, G.P. and P.A. Charette, 1975, The permittivity and attenuation
%   in polycrystalline and single-crystal ice Ih at 35 and 60 MHz, J.
%   Glaciol., 14(71), 293-303
%   
%   Stillman, D.E., J.A. MacGregor and R.E. Grimm, 2013, The role of acids
%   in electrical conduction through ice, J. Geophys. Res., 118
% 
%   Wolff, E.W., W.D. Miners, J.C. Moore and J.G. Paren, 1997, Factors
%   controlling the electrical conductivity of ice from the polar regions -
%   a summary, J. Phys. Chem. B, 101, 6090-6094.
%   
% Joe MacGregor, joemac@ig.utexas.edu
% Last updated: 02/25/15

function hfcp               = hfcondprop(varargin)

if (nargin && mod(nargin, 1))
    error('hfcondprop:nargin', 'Incorrect number of inputs (must be odd).')
end
if nargin
    if ~ischar(varargin{1})
        error('hfcondprop:modelstr', 'Input model must be a string.')
    elseif ~any(strcmp(varargin{1}, {'M07-std' 'M07-adj' 'M07-adj-NH4' 'M07-wt' 'W97' 'W97-adj'}))
        error('hfcondprop:modeltype', 'Input model must be either ''M07-std'', ''M07-adj'', ''M07-adj-NH4'', ''M07-wt'', ''W97'', ''W97-adj''.')
    else
        model               = varargin{1};
    end
    if (nargin > 1)
        for ii = 2:2:nargin
            if (~ischar(varargin{ii}) || ~any(strcmp(varargin{ii}, {'do_premelt' 'do_eutectic' 'activ_energy_pure' 'activ_energy_H' 'activ_energy_Cl' 'activ_energy_NH4' 'conduct_pure' 'molar_conduct_H' 'molar_conduct_Cl' 'molar_conduct_NH4' 'activ_energy_pure_uncert' 'activ_energy_H_uncert' ...
                                                                    'activ_energy_Cl_uncert' 'activ_energy_NH4_uncert' 'conduct_pure_uncert' 'molar_conduct_H_uncert' 'molar_conduct_Cl_uncert' 'molar_conduct_NH4_uncert' 'fact_premelt' 'temp_ref' 'temp_premelt' 'temp_eutectic_Cl'})))
                error('hfcondprop:paramname', 'Parameter to be adjusted is not a valid HFCP field.')
            end
            if (~(islogical(varargin{ii + 1}) || isnumeric(varargin{ii + 1})) || (length(varargin{ii + 1}) ~= 1))
                error('hfcondprop:paramtype', 'Assigned value for parameter to be adjusted is not a valid (must be a logical or numeric scalar).')
            end
        end
    end
else
    model                   = 'M07-std';
end
if (nargout > 1)
    error('hfcondprop:nargout', 'Incorrect number of outputs (must be 1).')
end

hfcp                        = struct;

% reported values at reference temperature (1st row) and their uncertainties (2nd row)
% activation energies: eV; conductivities: 10^-6 S/m; molar conductivities: S/m/(mol/L)
switch model
    case {'M07-std' 'M07-adj' 'M07-adj-NH4' 'M07-wt'}
        switch model 
            case {'M07-std' 'M07-wt'}
                activ_energy_pure ...
                            = [0.61 0.51 0.585 0.51; ...
                               NaN  0.01 0.024 0.01];
                conduct_pure= [4.5 9.2 6.0; ...
                               NaN 0.2 0.1];
            case {'M07-adj' 'M07-adj-NH4'}
                activ_energy_pure ...
                            = [0.51;
                               0.01];
                conduct_pure= [9.2;
                               0.2];
        end
        activ_energy_H      = [0.195 0.26 0.20 0.16; ...
                               NaN   0.03 0.01 NaN];
        activ_energy_Cl     = [0.19; ...
                               0.02];
        molar_conduct_H     = [2.57 3.5 3.66 3.5 2.3 3.0 3.1 3.8; ...
                               0.09 0.3 0.16 0.5 1.0 NaN 0.7 NaN];
        molar_conduct_Cl    = [0.43; ...
                               0.01];
    case {'W97' 'W97-adj'}
        activ_energy_pure   = [0.58; ...
                               0];
        activ_energy_H      = [0.21; ...
                               0];
        activ_energy_Cl     = [0.23; ...
                               0];
        conduct_pure        = [9; ...
                               0];
        molar_conduct_H     = [4;
                               0];
        molar_conduct_Cl    = [0.55; ...
                               0];
end
switch model
    case {'M07-std' 'M07-adj' 'M07-wt'}
        activ_energy_NH4    = [0; ...
                               0];
        molar_conduct_NH4   = [0; ...
                               0];
    case {'M07-adj-NH4' 'W97' 'W97-adj'}
        activ_energy_NH4    = [0.23; ...
                               0];
        switch model
            case 'M07-adj-NH4'
                molar_conduct_NH4 ...
                            = [0.8; ...
                               0];
            case {'W97' 'W97-adj'}
                molar_conduct_NH4 ...
                            = [1; ...
                               0];
        end
end
if strcmp(model, 'W97-adj')
    [conduct_pure, molar_conduct_H, molar_conduct_Cl, molar_conduct_NH4] ...
                            = deal((2.6 * conduct_pure), (2.6 * molar_conduct_H), (2.6 * molar_conduct_Cl), (2.6 * molar_conduct_NH4));
end

switch model
    
    case 'M07-std'
        
        [activ_energy_pure, activ_energy_pure_uncert, activ_energy_H, activ_energy_H_uncert, activ_energy_Cl, activ_energy_Cl_uncert, activ_energy_NH4, activ_energy_NH4_uncert, conduct_pure, conduct_pure_uncert, molar_conduct_H, molar_conduct_H_uncert, molar_conduct_Cl, molar_conduct_NH4, ...
         molar_conduct_NH4_uncert] ...
                            = deal(mean(activ_energy_pure(1, :)), std(activ_energy_pure(1, :)), mean(activ_energy_H(1, :)), std(activ_energy_H(1, :)), activ_energy_Cl(1), activ_energy_Cl(2), activ_energy_NH4(1), activ_energy_NH4(2), mean(conduct_pure(1, :)), std(conduct_pure(1, :)), ...
                                   mean(molar_conduct_H(1, :)), std(molar_conduct_H(1, :)), molar_conduct_Cl(1), molar_conduct_NH4(1), molar_conduct_NH4(2));
        molar_conduct_Cl_uncert ...
                            = (molar_conduct_H_uncert / molar_conduct_H) * molar_conduct_Cl(1); % assign molar_conduct_Cl same relative uncertainty as molar_conduct_H
        
    case 'M07-wt'
        
        activ_energy_pure(2, isnan(activ_energy_pure(2, :))) ...
                            = activ_energy_pure(1, isnan(activ_energy_pure(2, :))) .* mean(activ_energy_pure(2, ~isnan(activ_energy_pure(2, :))) ./ activ_energy_pure(1, ~isnan(activ_energy_pure(2, :))));
        activ_energy_pure_uncert ...
                            = sqrt(1 / sum(1 ./ (activ_energy_pure(2, :) .^ 2)));
        activ_energy_pure   = sum(activ_energy_pure(1, :) ./ (activ_energy_pure(2, :) .^ 2)) / sum(1 ./ (activ_energy_pure(2, :) .^ 2));
        activ_energy_H(2, isnan(activ_energy_H(2, :))) ...
                            = activ_energy_H(1, isnan(activ_energy_H(2, :))) .* mean(activ_energy_H(2, ~isnan(activ_energy_H(2, :))) ./ activ_energy_H(1, ~isnan(activ_energy_H(2, :))));
        activ_energy_H_uncert ...
                            = sqrt(1 / sum(1 ./ (activ_energy_H(2, :) .^ 2)));
        activ_energy_H      = sum(activ_energy_H(1, :) ./ (activ_energy_H(2, :) .^ 2)) / sum(1 ./ (activ_energy_H(2, :) .^ 2));
        activ_energy_Cl_uncert ...
                            = activ_energy_Cl(2);
        activ_energy_Cl     = activ_energy_Cl(1);
        activ_energy_NH4    = activ_energy_NH4(1);
        activ_energy_NH4    = activ_energy_NH4(2);
        conduct_pure(2, isnan(conduct_pure(2, :))) ...
                            = conduct_pure(1, isnan(conduct_pure(2, :))) .* mean(conduct_pure(2, ~isnan(conduct_pure(2, :))) ./ conduct_pure(1, ~isnan(conduct_pure(2, :))));
        conduct_pure_uncert = sqrt(1 / sum(1 ./ (conduct_pure(2, :) .^ 2)));
        conduct_pure        = sum(conduct_pure(1, :) ./ (conduct_pure(2, :) .^ 2)) / sum(1 ./ (conduct_pure(2, :) .^ 2));
        molar_conduct_H(2, isnan(molar_conduct_H(2, :))) ...
                            = molar_conduct_H(1, isnan(molar_conduct_H(2, :))) .* mean(molar_conduct_H(2, ~isnan(molar_conduct_H(2, :))) ./ molar_conduct_H(1, ~isnan(molar_conduct_H(2, :))));
        molar_conduct_H_uncert ...
                            = sqrt(1 / sum(1 ./ (molar_conduct_H(2, :) .^ 2)));
        molar_conduct_H  = sum(molar_conduct_H(1, :) ./ (molar_conduct_H(2, :) .^ 2)) / sum(1 ./ (molar_conduct_H(2, :) .^ 2));
        molar_conduct_Cl    = molar_conduct_Cl(1);
        molar_conduct_Cl_uncert ...
                            = (molar_conduct_H_uncert / molar_conduct_H) * molar_conduct_Cl;
        molar_conduct_NH4   = molar_conduct_NH4(1);
        molar_conduct_NH4_uncert ...
                            = molar_conduct_NH4(2);
        
    case {'M07-adj' 'M07-adj-NH4'}
        
        [activ_energy_pure, activ_energy_pure_uncert, activ_energy_H, activ_energy_H_uncert, activ_energy_Cl, activ_energy_Cl_uncert, activ_energy_NH4, activ_energy_NH4_uncert, conduct_pure, conduct_pure_uncert, molar_conduct_H, molar_conduct_H_uncert, molar_conduct_Cl, ...
         molar_conduct_Cl_uncert, molar_conduct_NH4, molar_conduct_NH4_uncert] ...
                            = deal(activ_energy_pure(1), activ_energy_pure(2), mean(activ_energy_H(1, :)), std(activ_energy_H(1, :)), activ_energy_Cl(1), activ_energy_Cl(2), activ_energy_NH4(1), activ_energy_NH4(2), conduct_pure(1), conduct_pure(2), mean(molar_conduct_H(1, :)), ...
                                   std(molar_conduct_H(1, :)), molar_conduct_Cl(1), molar_conduct_Cl(2), molar_conduct_NH4(1), molar_conduct_NH4(2));
        
    case {'W97' 'W97-adj'}
        
        [activ_energy_pure, activ_energy_pure_uncert, activ_energy_H, activ_energy_H_uncert, activ_energy_Cl, activ_energy_Cl_uncert, activ_energy_NH4, activ_energy_NH4_uncert, conduct_pure, conduct_pure_uncert, molar_conduct_H, molar_conduct_H_uncert, molar_conduct_Cl, ...
         molar_conduct_Cl_uncert, molar_conduct_NH4, molar_conduct_NH4_uncert] ...
                            = deal(activ_energy_pure(1), activ_energy_pure(2), activ_energy_H(1), activ_energy_H(2), activ_energy_Cl(1), activ_energy_Cl(2), activ_energy_NH4(1), activ_energy_NH4(2), conduct_pure(1), conduct_pure(2), molar_conduct_H(1), molar_conduct_H(2), ...
                                   molar_conduct_Cl(1), molar_conduct_Cl(2), molar_conduct_NH4(1), molar_conduct_NH4(2));
end

% load variables into hfcp structure and correct for significant figures
hfcp.activ_energy_H         = 1e-2 * round(1e2 * activ_energy_H);
hfcp.activ_energy_H_uncert  = 1e-2 * round(1e2 * activ_energy_H_uncert);
hfcp.activ_energy_Cl        = 1e-2 * round(1e2 * activ_energy_Cl);
hfcp.activ_energy_Cl_uncert = 1e-2 * round(1e2 * activ_energy_Cl_uncert);
hfcp.activ_energy_NH4       = 1e-2 * round(1e2 * activ_energy_NH4);
hfcp.activ_energy_NH4_uncert= 1e-2 * round(1e2 * activ_energy_NH4_uncert);
hfcp.conduct_pure           = 1e-7 * round(1e1 * conduct_pure);
hfcp.conduct_pure_uncert    = 1e-7 * round(1e1 * conduct_pure_uncert);
hfcp.molar_conduct_H        = 1e-1 * round(1e1 * molar_conduct_H);
hfcp.molar_conduct_H_uncert = 1e-1 * round(1e1 * molar_conduct_H_uncert);
hfcp.molar_conduct_Cl       = 1e-2 * round(1e2 * molar_conduct_Cl);
hfcp.molar_conduct_Cl_uncert= 1e-2 * round(1e2 * molar_conduct_Cl_uncert);
hfcp.molar_conduct_NH4      = 1e-2 * round(1e2 * molar_conduct_NH4);
hfcp.molar_conduct_NH4_uncert ...
                            = 1e-2 * round(1e2 * molar_conduct_NH4_uncert);

switch model
    case {'M07-std' 'M07-adj' 'M07-adj-NH4' 'M07-wt'}
        hfcp.activ_energy_pure ...
                            = 1e-2 * round(1e2 * activ_energy_pure);
        hfcp.activ_energy_pure_uncert ...
                            = 1e-2 * round(1e2 * activ_energy_pure_uncert);
    case {'W97' 'W97-adj'}
        [hfcp.activ_energy_pure, hfcp.activ_energy_pure_uncert] ...
                            = deal(activ_energy_pure, activ_energy_pure_uncert);
end

hfcp.fact_premelt           = 2; % premelting factor increasing conductivity at temp_premelt, assumed folowing MacGregor et al. [2007] (typically ignored)
switch model
    case {'M07-std' 'M07-adj' 'M07-adj-NH4' 'M07-wt'}
        hfcp.temp_ref       = 251; % reference temperature, K
    case {'W97' 'W97-adj'}
        hfcp.temp_ref       = 273.15 - 15;
end
hfcp.temp_premelt           = 263; % temperature at which pre-melting takes effect, assumed folowing MacGregor et al. [2007] (typically ignored)
hfcp.temp_eutectic_Cl       = 252; % eutectic temperature of chloride, -21 C, in K
switch model
    case {'M07-std' 'M07-adj' 'M07-wt'}
        hfcp.permitt_ice    = 3.2; % real part of complex relative permittivity of ice
    case {'M07-adj-NH4' 'W97' 'W97-adj'}
        hfcp.permitt_ice    = 3.15;
end
hfcp.do_premelt             = false; % true/false switch for raising conductivity above temp_premelt (ignored by default)
hfcp.do_eutectic_Cl         = false; % true/false switch for raising conductivity above temp_eutectic_Cl (ignored by default)
hfcp.model                  = model; % preserve model name

if (nargout > 1)
    hfcp.model              = [hfcp.model ' WITH ADDITIONAL ADJUSTMENTS'];
end

% adjust dielectric properties based on additional inputs during function call
if (nargin > 1)
    for ii = 2:2:nargin
        eval(['hfcp.' varargin{ii} ' = varargin{ii + 1};'])
    end
end

hfcp                        = orderfields(hfcp);
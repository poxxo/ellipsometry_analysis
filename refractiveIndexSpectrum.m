function n = refractiveIndexSpectrum(lambda, model, varargin)
%REFRACTIVEINDEXSPECTRUM Returns the refractive index of the specified
%material at desired wavelengths.
%   'lamda' is an array of wavelengths (in nm) at which the refractive 
%   index is returned.
%   'model' is either the name of a material, or a dispersion model.
%   If 'model' is a dispersion model, parameters of the model are specified
%   in additional arguments.
%
%%% Cauchy
%   Parameters: N0, N1, N2, ..., K0, K1, K2, ...
%   n = N0 + N1/lambda^2 + N2/lambda^4 + ...
%   k = K0 + K1/lambda^2 + K2/lambda^4 + ...
%
%%% Lorentz
%   Parameters: eOffset (dimensionless); and En (eV), Br (eV), Amp (dimensionless) for each oscillator.
%   Amp is approximately e2 at its peak value, Br is approximately the FWHM
%   of e2 and En is the resonant energy.
%   e = eOffset + (Amp*Br*En)/(En^2 - E^2 - i*Br*E) = e1 + 1j*e2
%
%%% Harmonic
%   Parameters: eOffset (dimensionless); and En (eV), Br (eV), Amp (dimensionless) for each oscillator.
%   Amp is approximately e2 at its peak value, Br is approximately the FWHM
%   of e2 and En is the resonant energy.
%   e = eOffset + 1/2*(Amp*Br)/(1/(En - E - i*Br/2) + 1/(En + E - i*Br/2)) = e1 + 1j*e2

hc_q = 6.62607004e-34*2.99792458e8/1.60217662e-19*1e9; % ~1240, eV-nm

n = zeros(1, length(lambda));
epsilon = zeros(1, length(lambda));
if strcmp(model, 'SiliconRegpro') || strcmp(model, 'SiliconPalik') || strcmp(model, 'SiO2Malitson') || strcmp(model, 'SiO2RodriguezDeMarcos')
    load(strcat('db/', strcat(model, '.mat')));
    n = interp1(lam, nData, lambda, 'spline') - 1j*interp1(lam, kData, lambda, 'spline');
elseif strcmp(model, 'Cauchy')
    numParams = (nargin - 2)/2;
    for i = 1:numParams
        n = n + varargin{i}./(lambda.^(2*(i - 1))) - 1j*varargin{numParams + i}./(lambda.^(2*(i - 1)));
    end
elseif strcmp(model, 'Lorentz')
    E = hc_q./lambda;
    numOsc = (nargin - 3)/3;
    for i = 1:numOsc
        epsilon = epsilon + (varargin{1 + 3*(i - 1) + 3}*varargin{1 + 3*(i - 1) + 2}*varargin{1 + 3*(i - 1) + 1})./(varargin{1 + 3*(i - 1) + 1}^2 - E.^2 - 1j*varargin{1 + 3*(i - 1) + 2}*E);
    end
    epsilon = epsilon + varargin{1};
    n = conj(sqrt(epsilon));
elseif strcmp(model, 'Harmonic')
    E = hc_q./lambda;
    numOsc = (nargin - 3)/3;
    for i = 1:numOsc
        epsilon = epsilon + 0.5*(varargin{1 + 3*(i - 1) + 3}*varargin{1 + 3*(i - 1) + 2}).*( 1./(varargin{1 + 3*(i - 1) + 1} - E - 0.5*1j*varargin{1 + 3*(i - 1) + 2}) + 1./(varargin{1 + 3*(i - 1) + 1} + E - 0.5*1j*varargin{1 + 3*(i - 1) + 2}) );
    end
    epsilon = epsilon + varargin{1};
    n = conj(sqrt(epsilon));
end

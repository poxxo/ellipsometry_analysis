function varargout = tlAnalysis( theta, lambda, n, L, varargin )
%TLANALYSIS Returns s- and p-polarization reflection coefficients or
%relevant parameters for a given layered structure.
%   theta: Array of angles (degrees)
%   lambda: Array of wavelengths (nm)
%   n: Refractive index (with n - jk convention) array of size
%   Nlayer-by-Nlambda. n can also be a vector of length Nlayer: in this
%   case, it is assumed that all layers have a wavelength-independent
%   refractive index.
%   L: Array of layer thicknesses (nm) of length Nlayer - 2
%   outputType (optional): String specifying the output type:
%   -- 'Reflectance' (default): s- and p-polarization reflectances,
%   returned as two Ntheta-by-Nlambda matrices (default)
%   -- 'ReflectionCoeff': s- and p-polarization reflection coefficients
%   -- 'ReflectionCoeffAll': s- and p-polarization reflection coefficients
%   at each interface, returned as two Ntheta-by-Nlambda-by-(Nlayer-1)
%   arrays.
%   -- 'PsiDelta': Psi and Delta (degrees), returned as two
%   Ntheta-by-Nlambda matrices, where tan(Psi)*exp(Delta) = r_p/r_s
%   -- 'TanCos': tan(Psi) and cos(Delta), returned as two Ntheta-by-Nlambda
%   matrices
%   -- 'ComplexPowerFlow': Complex powers at interfaces, normalized to the 
%   input real power, returned as two Ntheta-by-Nlambda-by-(Nlayer-1)
%   arrays.
%   -- 'PowerFlow': Real powers at interfaces, normalized to the
%   input real power, returned as two Ntheta-by-Nlambda-by-(Nlayer-1)
%   arrays.
%   -- 'ComplexPowerDep': Complex power deposited in each layer, normalized
%   to the  input real power, returned as two Ntheta-by-Nlambda-by-Nlayer
%   arrays. The matrix corresponding to the first layer is the complex
%   power reflected back.
%   -- 'PowerDep': Real power deposited in each layer, normalized to the
%   input real power, returned as two Ntheta-by-Nlambda-by-Nlayer arrays.
%   The matrix corresponding to the first layer is the real power
%   reflected back.
%   -- 'FieldCoeff': Returns the voltage coefficients A such that in each
%   layer, the voltage is of the form A(exp(g*L) + G*exp(-g*L)). For the
%   last layer, the voltage is of the form A*exp(-g*z). Two arrays of type
%   Ntheta-by-Nlambda-by-Nlayer are returned: one for each polarization.
%   -- 'TanCosCat': tan(Psi) and cos(Delta) concatenated along wavelength
%   direction, for fit purposes only.
%
%   Note: When theta consists of real elements, the first layer must have a
%   real refractive index for meaningful results.
%
%   Note: The choice of units of length for lambda and thickness does not
%   matter as long as they are consistent.

assert(nargin == 4 || nargin == 5, 'Expected 4 or 5 arguments');
if nargin == 4
    outputType = 'Reflectance';
else
    assert(     strcmpi(varargin{1}, 'Reflectance') ...
             || strcmpi(varargin{1}, 'ReflectionCoeff') ...
             || strcmpi(varargin{1}, 'ReflectionCoeffAll') ...
             || strcmpi(varargin{1}, 'PsiDelta') ...
             || strcmpi(varargin{1}, 'TanCos') ...
             || strcmpi(varargin{1}, 'TanCosCat') ...
             || strcmpi(varargin{1}, 'ComplexPowerFlow') ...
             || strcmpi(varargin{1}, 'PowerFlow') || strcmpi(varargin{1}, 'RealPowerFlow') ...
             || strcmpi(varargin{1}, 'ComplexPowerDep') ...
             || strcmpi(varargin{1}, 'PowerDep') || strcmpi(varargin{1}, 'RealPowerDep') ...
             || strcmpi(varargin{1}, 'FieldCoeff') ...
         , 'Invalid output type');
    outputType = varargin{1};
end

assert(isvector(theta), 'theta must be a vector');
assert(isvector(lambda), 'lambda must be a vector');
assert(isvector(L) || isempty(L), 'L must be a vector or an empty array');
assert((size(n, 2) == length(lambda) && size(n, 1) >= 2) || (isvector(n) && length(n) >= 2), ...
    'refInx has invalid dimensions: expected either Nlayer by Nlambda or a vector of length Nlayer (Nlayer >= 2)')

if size(lambda, 1) > 1
    lambda = lambda.';
end
if size(theta, 1) > 1
    theta = theta.';
end
if isvector(n) % rn is constant in wavelength
    if size(n, 2) > 1
        n = n.';
    end
    n = repmat(n, 1, length(lambda));
end

assert(length(L) + 2 == size(n, 1), 'n and L have incompatible dimensions');


lambda = reshape(lambda, 1, length(lambda));
theta = reshape(theta, 1, length(theta));

Ntheta = length(theta);
Nlambda = length(lambda);
Nlayer = size(n, 1); % Number of layers, including the first and last


% Cosines of the transmission angles
cosTrans = zeros(Ntheta, Nlambda, Nlayer);
cosTrans(:, :, 1) = transpose(cosd(theta))*ones(1, Nlambda);
%signumUp = @(x) sign(x) + (x == 0);
for i = 1:(Nlayer - 1)
    %cosTrans(:, :, i + 1) = signumUp(-imag(sqrt(1 - (transpose(sind(theta))*(n(1, :)./n(i + 1, :))).^2))) .* sqrt(1 - (transpose(sind(theta))*(n(1, :)./n(i + 1, :))).^2);
    cosTrans(:, :, i + 1) = sqrt(1 - (transpose(sind(theta))*(n(1, :)./n(i + 1, :))).^2);
end

% All layers are implicitly assumed to consist of non-magnetic materials (relative permeability = 1).
% If this is not the case, n should be replaced with the reciprocal of the wave impedance sqrt(mu_r/epsilon_r).
n = repmat(permute(n, [3, 2, 1]), Ntheta, 1, 1);

% Propagation constants
g = 1j*2*pi*n.*cosTrans./repmat(lambda, Ntheta, 1, Nlayer);

% Characteristic impedances
Zch_s = 1./(n.*cosTrans);
Zch_p = cosTrans./n;

% Interface impedances
Z_s = zeros(Ntheta, Nlambda, Nlayer - 1);
Z_p = zeros(Ntheta, Nlambda, Nlayer - 1);
Z_s(:, :, Nlayer - 1) = Zch_s(:, :, Nlayer);
Z_p(:, :, Nlayer - 1) = Zch_p(:, :, Nlayer);
if Nlayer > 2
    for i = linspace(Nlayer - 2, 1, Nlayer - 2)
        Z_s(:, :, i) = (Z_s(:, :, i+1) + Zch_s(:, :, i+1).*tanh(g(:, :, i+1)*L(i))) ./ (Zch_s(:, :, i+1) + Z_s(:, :, i+1).*tanh(g(:, :, i+1)*L(i))) .* Zch_s(:, :, i+1);
        Z_p(:, :, i) = (Z_p(:, :, i+1) + Zch_p(:, :, i+1).*tanh(g(:, :, i+1)*L(i))) ./ (Zch_p(:, :, i+1) + Z_p(:, :, i+1).*tanh(g(:, :, i+1)*L(i))) .* Zch_p(:, :, i+1);
    end
end

% Reflection coefficients
G_s = zeros(Ntheta, Nlambda, Nlayer - 1);
G_p = zeros(Ntheta, Nlambda, Nlayer - 1);
for i = 1:(Nlayer - 1)
    G_s(:, :, i) = (Z_s(:, :, i) - Zch_s(:, :, i))./(Z_s(:, :, i) + Zch_s(:, :, i));
    G_p(:, :, i) = (Z_p(:, :, i) - Zch_p(:, :, i))./(Z_p(:, :, i) + Zch_p(:, :, i));
end
if strcmpi(outputType, 'ReflectionCoeff')
    varargout{1} = G_s(:, :, 1);
    varargout{2} = G_p(:, :, 1);
    return;
end
if strcmpi(outputType, 'ReflectionCoeffAll')
    varargout{1} = G_s;
    varargout{2} = G_p;
    return;
end

% G0_s = G_s(:, :, 1); % Reflection coefficients at the first boundary
% G0_p = G_p(:, :, 1);
% R0_s = abs(G0_s).^2 + 2*imag(G0_s).*imag(Zch_s(:, :, 1))./real(Zch_s(:, :, 1)); % Reflectances
% R0_p = abs(G0_p).^2 + 2*imag(G0_p).*imag(Zch_p(:, :, 1))./real(Zch_p(:, :, 1));

if strcmpi(outputType, 'Reflectance')
    % The definition of reflectance requires care when the medium from
    % which the wave is incident is lossy (i.e. has a complex
    % characteristic impedance)
    % Here we take the definition to be "the real component of the
    % reflected complex power, divided by the real component of the
    % incident complex power (the power if the input reflection coefficient
    % were zero)"
    varargout{1} = abs(G_s(:, :, 1)).^2 + 2*imag(G_s(:, :, 1)).*imag(Zch_s(:, :, 1))./real(Zch_s(:, :, 1));
    varargout{2} = abs(G_p(:, :, 1)).^2 + 2*imag(G_p(:, :, 1)).*imag(Zch_p(:, :, 1))./real(Zch_p(:, :, 1));
    return;
end


if strcmpi(outputType, 'PsiDelta')
    r_s = G_s(:, :, 1);
    r_p = -G_p(:, :, 1);
    % Handedness flips upon reflection. The negative sign is required by the field vector orientations assumed by above equations.
    % This is a matter of convention. With the other convention, r_s = 1 and r_p = 1 gives RHEP -> LHEP and vice versa.
    varargout{1} = atand(abs(r_p./r_s));        % psi
    varargout{2} = 180.0/pi*angle(r_p./r_s);    % delta
    return;
end
if strcmpi(outputType, 'TanCos')
    r_s = G_s(:, :, 1);
    r_p = -G_p(:, :, 1);
    varargout{1} = abs(r_p./r_s);           % tan(psi)
    varargout{2} = cos(angle(r_p./r_s));    % cos(delta)
    return;
end
if strcmpi(outputType, 'TanCosCat')
    r_s = G_s(:, :, 1);
    r_p = -G_p(:, :, 1);
    varargout{1} = cat(2, abs(r_p./r_s), cos(angle(r_p./r_s)));
    return;
end

%%%%

% Field expression coefficients
A_s = zeros(Ntheta, Nlambda, Nlayer);
A_p = zeros(Ntheta, Nlambda, Nlayer);
A_s(:, :, 1) = sqrt(abs(2./real(1./Zch_s(:, :, 1)))); % Normalized such that the incident power is unity
A_p(:, :, 1) = sqrt(abs(2./real(1./Zch_p(:, :, 1))));
if Nlayer > 2
    for i = 2:(Nlayer - 1)
        A_s(:, :, i) = A_s(:, :, i - 1) .* (1 + G_s(:, :, i - 1)) ./ (exp(g(:, :, i)*L(i-1)) + G_s(:, :, i).*exp(-g(:, :, i)*L(i-1)));
        A_p(:, :, i) = A_p(:, :, i - 1) .* (1 + G_p(:, :, i - 1)) ./ (exp(g(:, :, i)*L(i-1)) + G_p(:, :, i).*exp(-g(:, :, i)*L(i-1)));
    end
end
A_s(:, :, Nlayer) = A_s(:, :, Nlayer - 1) .* (1 + G_s(:, :, Nlayer - 1));
A_p(:, :, Nlayer) = A_p(:, :, Nlayer - 1) .* (1 + G_p(:, :, Nlayer - 1));

if strcmpi(outputType, 'FieldCoeff')
    varargout{1} = A_s;
    varargout{2} = A_p;
    return;
end

% Complex power flows at interfaces
S_s = abs((1 + G_s).*A_s(:, :, 1:(Nlayer - 1))).^2./(2*conj(Z_s));
S_p = abs((1 + G_p).*A_p(:, :, 1:(Nlayer - 1))).^2./(2*conj(Z_p));
% S_s = abs(A_s(:, :, 1:(Nlayer - 1))).^2.*(1 + G_s).*conj(1 - G_s)./(2*conj(Zch_s(:, :, 1:(Nlayer - 1))));
% S_p = abs(A_p(:, :, 1:(Nlayer - 1))).^2.*(1 + G_p).*conj(1 - G_p)./(2*conj(Zch_p(:, :, 1:(Nlayer - 1))));
if strcmpi(outputType, 'ComplexPowerFlow')
    varargout{1} = S_s;
    varargout{2} = S_p;
    return;
end
if strcmpi(outputType, 'PowerFlow') || strcmpi(outputType, 'RealPowerFlow')
    varargout{1} = real(S_s);
    varargout{2} = real(S_p);
    return;
end

% Complex powers deposited in layers
% P_s = cat(3, 1 - S_s(:, :, 1), S_s - cat(3, zeros(Ntheta, Nlambda), S_s(:, :, 2:end)));
% P_p = cat(3, 1 - S_p(:, :, 1), S_p - cat(3, zeros(Ntheta, Nlambda), S_p(:, :, 2:end)));
P_s = cat(3, 1 - S_s(:, :, 1), S_s - cat(3, S_s(:, :, 2:end), zeros(Ntheta, Nlambda)));
P_p = cat(3, 1 - S_p(:, :, 1), S_p - cat(3, S_p(:, :, 2:end), zeros(Ntheta, Nlambda)));
if strcmpi(outputType, 'ComplexPowerDep')
    varargout{1} = P_s;
    varargout{2} = P_p;
    return;
end
if strcmpi(outputType, 'PowerDep') || strcmpi(outputType, 'RealPowerDep')
    varargout{1} = real(P_s);
    varargout{2} = real(P_p);
    return;
end

error('Unexpected output type');

end


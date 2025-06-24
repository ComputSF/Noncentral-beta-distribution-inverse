# Noncentral-beta-distribution-inverse
Inversion of the noncentral beta distribution function

# 1. Identification of the noncentrality parameter


`xinv` is a MATLAB function designed to compute the noncentrality parameter of the Noncentral Beta Distribution. 

## Function Signature

```matlab
[x,ier] = xinv(z, y, p, q)
```

## Inputs

- **`z`**: Target cumulative probability (should be in the range [0, 1])
- **`y`**: Quantile or value at which the noncentral beta CDF is evaluated (should be in the range [0, 1])
- **`p`**: Shape parameter \( p \) (must be positive)
- **`q`**: Shape parameter \( q \) (must be positive)

## Outputs

- **`x`**: Computed value of the noncentrality parameter (lambda)
- **`ierr`**: Error flag indicating the success or issues in computation:
  - `0`: Computation successful
  - `1`: Some loss of accuracy is expected
  - `2`: The z-value is not acceptable for the given y, p, q values: z > I_y(p,q). A negative value for x is returned.
  - `3`: The initial approximation to x is greater than 1000 (upper limit for the computation of Bpqxy). The initial approximation to x is returned.

## Description

The algorithm combines sophisticated asymptotic approximations for generating accurate initial values with the highly efficient Halley's method for iterative refinement. 

## Usage

Here is a basic example of how to use the `xinv` function:

```matlab
% Given parameters
z = 0.0043; % Cumulative probability
y = 0.12;   % Input value (should be in the range [0, 1])
p = 3.0;    % Shape parameter p (must be positive)
q = 4.0;    % Shape parameter q (must be positive)

% Compute the noncentrality parameter x
[x, ier] = xinv(z, y, p, q);

% Display the results
if ier == 0
    fprintf('Computation successful.\n');
    fprintf('The noncentrality parameter x is: %f\n', x);
    
    % --- Verification ---
    z_computed = Bpqxy(x, y, p, q);
    fprintf('Verification: Bpqxy(x, y, p, q) = %f\n', z_computed);
    fprintf('Original z = %f\n', z);
    fprintf('Relative error = %e\n', abs(z_computed - z) / z);
    
else
    fprintf('Computation finished with an error flag.\n');
    fprintf('ier = %d\n', ier);
    fprintf('See documentation for the meaning of the error flag.\n');
    if ier == 3
        fprintf('The returned value is an initial approximation for x.\n');
        fprintf('x = %f\n', x);
    end
end


```


# 2. Calculation of the quantiles of the distribution


`yinv` is a MATLAB function designed to compute the quantiles of the Noncentral Beta Distribution. 

## Function Signature

```matlab
[y,ier]=yinv(z, x, p, q)
```

## Inputs

- **`z`**: Target cumulative probability (should be in the range [0, 1])
- **`x`**: Noncentrality parameter (lambda, should be positive)
- **`p`**: Shape parameter \( p \) (must be positive)
- **`q`**: Shape parameter \( q \) (must be positive)

## Outputs

- **`y`**: Computed value of the quantile
- **`ierr`**: Error flag indicating the success or issues in computation:
  - `0`: Computation successful
  - `1`: Some loss of accuracy is expected

## Description

The algorithm combines asymptotic approximations for generating accurate initial values with the highly efficient iterative refinement. 

## Usage

Here is a basic example of how to use the `yinv` function:

```matlab
% Given parameters
z = 0.0043; % Cumulative probability
x = 4.5;    % Noncentrality
p = 3.0;    % Shape parameter p (must be positive)
q = 4.0;    % Shape parameter q (must be positive)

% Compute the quantile y:
[y, ier] = yinv(z, x, p, q);

fprintf('The quantile y is: %f\n', y);

% --- Verification ---
z_computed = Bpqxy(x, y, p, q);
fprintf('Verification: Bpqxy(x, y, p, q) = %f\n', z_computed);
fprintf('Original z = %f\n', z);
fprintf('Relative error = %e\n', abs(z_computed - z) / z);


```

## Authors

V. Egorova, A. Gil, J. Segura and N. M. Temme

## Contact

For questions or comments, please contact Amparo Gil at [amparo.gil@unican.es].


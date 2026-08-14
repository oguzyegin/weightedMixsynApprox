# Weighted Mixed-Sensitivity Controller Approximation

MATLAB research code for approximating a high-/infinite-dimensional $H_\infty$ controller by a fixed-order rational controller while preserving closed-loop stability and minimizing a weighted mixed-sensitivity objective.

The repository combines an infinite-dimensional time-delay $H_\infty$ design with several fixed-order approximation and tuning approaches so their achieved mixed-sensitivity levels and computational costs can be compared on the same example.

## Problem

For a plant $P$, controller $C$, sensitivity

$$
S = (1 + PC)^{-1},
$$

and complementary sensitivity

$$
T = PC(1 + PC)^{-1},
$$

the performance measure used throughout the example is

$$
\gamma(C) =
\sup_{\omega}
\sqrt{
\left|W_1(j\omega)S(j\omega)\right|^2
+
\left|W_2(j\omega)T(j\omega)\right|^2
}.
$$

The goal is to obtain a low-order controller $C_m$ that stabilizes the plant and keeps this quantity as small as possible relative to the high-/infinite-dimensional reference controller.

## Main workflow

`main_example.m` runs the complete benchmark:

1. Defines the time-delay plant and the weights `W1` and `W2`.
2. Converts the numerator and denominator quasi-polynomials to a YALTA-compatible representation with `nPdP2YALTA`.
3. Detects unstable poles and constructs the inner/outer factorization required by the time-delay $H_\infty$ design.
4. Computes an infinite-dimensional reference controller `Copt` using `TO_hinfsyn`.
5. Samples the plant and controller frequency responses.
6. Generates fixed-order initial controllers using weighted state-space estimation and vector fitting.
7. Applies several fixed-order approximation/tuning methods.
8. Evaluates every controller using the same mixed-sensitivity performance measure.
9. Reports the achieved peak performance and computation time.

The supplied example uses a desired controller order of 8.

## Included approximation methods

### Weighted state-space estimation

`Sequential Convex Optimization/weightedSsest.m`

Fits an order-m state-space model to the reference-controller FRD. The fitting weight used in the main example is based on

$$
\left|W_x P S_{\mathrm{opt}}^2\right|,
\qquad
W_x = \sqrt{|W_1|^2 + |W_2|^2}.
$$

This controller is also used as an initializer for other methods.

### Vector fitting

`vf_fit_order8.m` and `vf_tune_hinf.m`

Constructs a stable fixed-order rational approximation and optionally tunes it directly against the delayed-plant mixed-sensitivity objective.

### `systune`

`callSystune.m`

Uses a finite-dimensional approximation of the plant and initializes a tunable transfer function from a supplied controller. The tuning objective is the weighted mixed-sensitivity channel, with an additional controller-stability requirement.

### TF-IRKA

`tfirkaMixedSensitivity.m`

Optimizes TF-IRKA interpolation points to construct an exactly order-$m$ controller. Candidate controllers are evaluated using the mixed-sensitivity objective, while a finite-dimensional plant approximation `Pa` is used for the closed-loop stability test.

Automatic initialization can use unstable plant poles and sensitivity peaks; user-specified interpolation points are also supported.

### Sequential convex optimization

`Sequential Convex Optimization/structuredMixedSensitivitySingleCluster.m`

The proposed fixed-order approximation method. It uses the reference-oriented residual

$$
E_{\mathrm{app}} =
W_x P(C_{\mathrm{opt}}-C)S S_{\mathrm{opt}}
$$

and a sequential SOCP procedure guided at active approximation-error peaks by the structure

$$
T_z(s) = (s-z)^r.
$$

The method supports a prescribed controller order and relative order. YALMIP and SDPT3 are used for the convex subproblems, while a rational plant approximation is used for closed-loop stability checks.

### TDS-CONTROL / HANSO

`designTDSMixedSensitivity8.m`

Directly optimizes an order-8, relative-degree-one controller on a dense frequency grid using the exact analytical delayed frequency response. TDS-CONTROL is used to test characteristic roots of the exact delayed closed loop, and unstable trial controllers are penalized.

## Reference time-delay controller

The directory `Hinf Optimal with YALTA v3/` contains the routines used to construct the high-/infinite-dimensional $H_\infty$ reference controller.

`nPdP2YALTA.m` converts delayed transfer-function numerator and denominator expressions into quasi-polynomial coefficient/delay data suitable for YALTA.

`TO_hinfsyn.m` interfaces with the supplied YALTA-based routines to compute the reference controller.

## Performance evaluation

`perfLevel.m` evaluates

$$
\sqrt{|W_1S|^2 + |W_2T|^2}
$$

over a specified frequency grid.

At the end of `main_example.m`, the controllers are compared using a MATLAB table containing:

- method name,
- maximum sampled mixed-sensitivity level,
- computation time.

The benchmark currently includes weighted `ssest`, vector fitting, TF-IRKA, `systune`, TDS-CONTROL/HANSO, and the sequential convex method.

## Repository structure

```text
main_example.m
nPdP2YALTA.m
TO_hinfsyn.m
perfLevel.m
callSystune.m
tfirkaMixedSensitivity.m
vf_fit_order8.m
vf_tune_hinf.m
designTDSMixedSensitivity8.m

Sequential Convex Optimization/
    structuredMixedSensitivitySingleCluster.m
    weightedSsest.m

Hinf Optimal with YALTA v3/
    YALTA-based time-delay H-infinity design routines

```

## Requirements

The code is written for MATLAB. Depending on the method being run, it uses:

- Control System Toolbox
- System Identification Toolbox
- Robust Control Toolbox (`systune`)
- Symbolic Math Toolbox (`nPdP2YALTA`)
- YALTA, included under `Hinf Optimal with YALTA v3/`
- TDS-CONTROL
- YALMIP and SDPT3 for the sequential convex method
- HANSO for `designTDSMixedSensitivity8`

The optional `particleswarm` mode of `tfirkaMixedSensitivity` requires Global Optimization Toolbox.

Some external packages may need to be added to the MATLAB path before their corresponding method is used.

## Running the example

Set the MATLAB working directory to the repository root and run:

```matlab
main_example
```

The script adds the supplied YALTA, sequential-convex-optimization, and TDS-CONTROL directories to the MATLAB path where needed.

For the sequential convex method, make sure the YALMIP directory and an SDPT3 installation are available on the MATLAB path. For the TDS-CONTROL optimization, TDS-CONTROL and HANSO packages must also be available.

## Notes

This repository is research code intended for experimentation with fixed-order approximation of mixed-sensitivity $H_\infty$ controllers for time-delay systems. Several methods use different finite-dimensional or frequency-domain representations for optimization and stability verification; the final performance comparison in `main_example.m` is evaluated against the original delayed plant on a dense frequency grid.

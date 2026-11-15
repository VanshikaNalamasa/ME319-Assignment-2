# ME319 Assignment 2: Predictive Tracking and Turret Elevation Control

**Course:** ME319 – Control Systems  
**Institution:** Indian Institute of Technology Bombay  
**Department:** Mechanical Engineering  
**Author:** Vanshika Nalamasa  
**Roll Number:** 24B2106  
**Semester:** Autumn 2026

---

## 1. Project Overview

This project implements a predictive-feedback control system for the elevation axis of a laser turret tasked with tracking a maneuvering aerial target.

The system uses noisy and intermittently unavailable radar measurements to estimate the target's position, velocity, and acceleration. The estimated target state is extrapolated over a prediction horizon of 0.5 seconds to generate a future elevation reference for the turret.

An Extended Kalman Filter (EKF), a Cartesian constant-acceleration prediction model, and a proportional controller are integrated to achieve predictive tracking. A separate uncertainty-aware firing mechanism determines when the estimated target position satisfies the engagement conditions.

The system is evaluated through three simulation scenarios involving non-maneuvering motion, target maneuvering, and radar measurement blackout.

## 2. System Architecture

The implemented system consists of five major stages:

1. **Radar measurement generation:** Generates noisy range and elevation measurements at 20 Hz, including measurement loss during the specified blackout interval.
2. **Target-state estimation:** Uses an Extended Kalman Filter to estimate the target's Cartesian position, velocity, and acceleration.
3. **Target prediction:** Extrapolates the estimated state over the 0.5-second prediction horizon.
4. **Turret control:** Uses a proportional controller to generate the motor voltage command required to track the predicted elevation.
5. **Firing logic:** Evaluates the estimated cross-beam distance and its uncertainty to determine whether a firing pulse should be issued.

The firing logic is independent of the voltage-control feedback loop.

## 3. Methodology

### 3.1 Target-State Estimation

The estimator uses the Cartesian state vector:

\[
\mathbf{s} =
[x,\ y,\ v_x,\ v_y,\ a_x,\ a_y]^T
\]

A constant-acceleration state transition model is used, with white-jerk process noise.

The EKF uses the nonlinear radar measurement model:

\[
h(\mathbf{s}) =
\begin{bmatrix}
\sqrt{x^2+y^2}\\
\operatorname{atan2}(y,x)
\end{bmatrix}
\]

The filter performs prediction at every controller step and measurement correction only when a new valid radar measurement is available.

During a radar blackout, the filter continues open-loop state propagation while its covariance increases.

### 3.2 Target Prediction

The estimated target state is propagated forward by:

\[
T_p = 0.5\text{ s}
\]

The predicted Cartesian position is converted into an elevation reference:

\[
\theta_{\mathrm{ref}} =
\operatorname{atan2}(y_p,x_p)
\]

The reference is clamped to the turret's permitted angular range.

### 3.3 Turret Plant and Controller

The turret is modelled as a third-order electromechanical system, incorporating motor electrical dynamics and mechanical rotational dynamics.

The equivalent inertia is:

\[
J_{\mathrm{eq}} =
J_{\mathrm{turret}}+N^2J_m
=3200\text{ kg m}^2
\]

The proportional controller is:

\[
V_m=K_p(\theta_{\mathrm{ref}}-\theta_{\mathrm{turret}})
\]

The controller gain is selected analytically using the prediction horizon and the plant velocity constant:

\[
K_p=\frac{1}{T_pK_v}
\]

The selected gain is approximately:

\[
K_p=2003\text{ V/rad}
\]

Stability is assessed using the Routh-Hurwitz criterion, closed-loop poles, gain margin, and phase margin.

### 3.4 Firing Logic

The firing mechanism uses the estimated present target position and turret orientation.

It evaluates the perpendicular distance between the target estimate and the turret beam, together with the estimated cross-beam uncertainty.

A firing pulse is issued only when:

- The target is in front of the turret.
- The uncertainty-aware engagement condition is satisfied.
- The minimum one-second pulse separation is maintained.
- The firing command transitions from low to high.

Each firing command lasts exactly one controller time step.

## 4. Simulation Scenarios

Three scenarios are implemented.

| Scenario | Target motion | Radar availability |
|---|---|---|
| A | Non-maneuvering target | Normal |
| B | Maneuvering target | Normal |
| C | Same maneuver as B | 0.5-second blackout |

The target trajectories are engineering assumptions implemented in the target-motion model. They are used for self-validation and should not be interpreted as instructor-provided evaluation trajectories.

## 5. Repository Structure

```text
ME319-Assignment-2/
│
├── README.md
│
├── params.m
├── turret_plant.m
├── target_truth.m
├── target_motion_generator.m
├── target_state_estimator.m
├── target_predictor.m
├── turret_controller.m
│
├── tests/
│   ├── test_stage2m.m
│   ├── test_turret_controller.m
│   └── [other estimator and generator tests]
│
├── figures/
│   ├── tmg_fig1_truth.png
│   ├── tmg_fig2_measurements.png
│   ├── fig1_elevation_tracking.png
│   ├── fig2_pointing_miss.png
│   ├── fig3_control_effort.png
│   ├── fig4_blackout_zoom.png
│   └── fig5_linear_loop.png
│
└── report/
    ├── ME319_Assignment2_Report.pdf
    └── main.tex
```

*The filenames and directory structure above should match the actual files uploaded to the repository.*

## 6. Software Requirements

- MATLAB
- MATLAB Control System Toolbox, if required by the linear loop analysis or test scripts.
- All project MATLAB source files and helper functions.

The exact toolbox requirements depend on which optional analysis scripts are executed.

## 7. Running the Project

1. Clone or download this repository.
2. Open MATLAB.
3. Set the MATLAB current folder to the repository root.
4. Ensure that the required source files and helper functions are on the MATLAB path.
5. Run the appropriate simulation or test script.

For example:

```matlab
clear;
clear functions;

test_turret_controller
```

Run the target-generator, estimator, and parameter/plant test scripts separately using their corresponding filenames.

**Note:** Use the actual test-script filenames present in the repository. The test scripts should be run from a clean MATLAB session to avoid dependence on previously loaded variables or functions.

## 8. Validation Results

The following results were obtained from the project's own MATLAB validation and closed-loop simulations.

### 8.1 Automated Test Results

| Component | Reported test result |
|---|---:|
| Parameter and plant validation | Checks passed |
| Target-motion generator | 63/63 |
| State estimator | 47/47 |
| Turret controller | 78/78 structural/design checks |
| Controller performance flags | 6/6 |

These results represent self-validation using the project's test scripts.

### 8.2 Closed-Loop Simulation Results

| Metric | Scenario A | Scenario B | Scenario C |
|---|---:|---:|---:|
| Pulses issued | 8 | 7 | 7 |
| Successful hits | 8 | 5 | 5 |
| First successful hit (s) | 2.700 | 4.200 | 4.130 |
| Miss RMS (m) | 0.77 | 1.21 | 1.19 |
| Maximum miss (m) | 1.27 | 2.62 | 2.54 |

The non-maneuvering scenario achieved successful engagement for all issued pulses. The maneuvering scenarios showed reduced pointing accuracy, with two missed engagements in each case.

The blackout scenario demonstrated continued state propagation during measurement loss, accompanied by increased estimation uncertainty.

## 9. Limitations

The implementation has the following limitations:

- The target estimator and predictor use a constant-acceleration model, which cannot perfectly represent changing target acceleration.
- Prediction and tracking accuracy degrade during maneuvering.
- Radar blackout causes estimation uncertainty to increase until valid measurements resume.
- The turret plant assumes an ideal linear drive.
- Voltage and current saturation are not modelled.
- Motor and amplifier hardware ratings are not verified.
- The target trajectories are self-defined engineering assumptions.
- The instructor's final evaluation has not yet been performed.

The peak simulated voltage and current are 645 V and 101 A, respectively. These are ideal-model demands and should not be interpreted as validated hardware operating limits.

## 10. Report

The detailed technical report contains:

- System architecture and mathematical modelling.
- EKF formulation and target prediction.
- Turret transfer-function derivation.
- Controller gain selection and stability analysis.
- Firing logic and uncertainty handling.
- Simulation results and scenario comparisons.
- Radar blackout analysis.
- Control effort and limitations.
- AI and computational tools disclosure.

The final PDF report and its LaTeX source are included in the `report/` directory.

## 11. Evaluation Status

The reported numerical

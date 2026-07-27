# OOK Biological Attenuation Compatibility Design

**Date:** 2026-07-27  
**Status:** Approved design  
**Reference implementation:** `krakenFortran`  
**Scope:** End-to-end support across the C++ core, ENV, JSON, public C++ API, Python, MATLAB, regression fixtures, Fortran differential checks, and TL acceptance.

## 1. Objective

Add the `B` biological attenuation option to OOK with the same valid-input numerical semantics as `krakenFortran`.

The implementation must:

- preserve the Fortran biological attenuation formula and operation order;
- preserve closed depth intervals and ordered overlap accumulation;
- transport biological-layer configuration through both `env_in_out.cpp` and `json_in_out.hpp`;
- expose the configuration through C++, Python, and MATLAB;
- reject malformed or physically undefined input before mutating active parameters;
- keep the current weak-loss complex-sound-speed conversion and SSP interpolation behavior;
- verify the implementation against `krakenFortran` at formula, modal-wavenumber, and transmission-loss levels.

The feature is incomplete if the core formula works but either the ENV or JSON transport path is missing.

## 2. Approved Decisions

The following decisions were explicitly approved:

1. Use an end-to-end implementation rather than a core-only patch.
2. Preserve Fortran behavior for valid input while making invalid input fail fast.
3. Keep biological attenuation as one global attenuation configuration shared by all profiles.
4. Reject multi-profile inputs whose attenuation configurations differ.
5. Extend the existing `Atten_Mode` object rather than adding parallel state.
6. Accept the C++ ABI change and rebuild all native OOK artifacts together.
7. Do not change the weak-loss complex-sound-speed formula.
8. Do not insert biological-layer boundary nodes automatically.
9. Treat `env_in_out.cpp` and `json_in_out.hpp` as mandatory, separately tested deliverables.

## 3. Existing Compatibility Gap

`krakenFortran` supports three mutually exclusive added volume-attenuation options in the fourth Top Option character:

- `T`: Thorp;
- `F`: Francois-Garrison;
- `B`: biological attenuation.

OOK currently supports only `None`, `Thorpe`, and `FrancGarr`. An ENV fourth character of `B` falls through to `None`, and the following biological-layer records are not consumed. The next line is then misinterpreted as either top-halfspace data or an internal-medium header.

OOK also has no JSON model, public parameter type, Python binding, or MATLAB helper capable of representing the biological layers.

## 4. Architecture

### 4.1 Single ownership

`OOK_parameters::AttenUnit` remains the sole owner of attenuation configuration.

No Fortran-style module global such as `bio[]` or `NBioLayers` will be introduced. This preserves OOK's object-based parameter flow, thread safety, reentrancy, and isolation between solver instances.

### 4.2 Public data model

Add the following public type near the existing attenuation types in `include/OpenOceanKrakenParams.h`:

```cpp
inline constexpr std::size_t MaxBioLayers = 200;

struct BiologicalAttenuationLayer
{
    double Z1 = 0.0;  // m
    double Z2 = 0.0;  // m
    double f0 = 0.0;  // Hz
    double Q  = 0.0;  // dimensionless
    double a0 = 0.0;  // dB/km
};
```

Append the new enum value rather than inserting it between existing values:

```cpp
enum class OceanAbsorptionModel
{
    None,
    Thorpe,
    FrancGarr,
    Biological
};
```

Extend the existing configuration atomically:

```cpp
struct Atten_Mode
{
    AttenuationUnit attnUnit =
        AttenuationUnit::MODE_W_db_per_lambda;

    OceanAbsorptionModel absModel =
        OceanAbsorptionModel::None;

    std::vector<BiologicalAttenuationLayer> biologicalLayers;
};
```

The layer count is always derived from `biologicalLayers.size()`. No separate `NBioLayers` field is stored.

### 4.3 Global multi-profile semantics

The biological-layer list belongs to the global `Atten_Mode`, not to an individual SSP medium or `Range_Independent_Area`.

All profiles in a multi-profile input must have identical:

- attenuation units;
- absorption model;
- biological-layer count;
- biological-layer values in the same order.

The parser must reject differences. It must not silently retain only the first profile's configuration.

### 4.4 Function and ownership flow

The existing flow remains:

```text
OOK_parameters::AttenUnit
    → input_SSP
    → UpdateSSPLoss
    → CRCI
    → parseAttenuation
    → addOceanAbsorption
    → complex cp/cs
```

`Atten_Mode` must be passed as `const Atten_Mode&` through the loss-update and attenuation functions so the vector is not copied.

The existing `Interface::set_AttenUnit()` remains the main public C++ and Python-native setter. No second independently mutable biological configuration is introduced.

## 5. ENV Transport

### 5.1 Canonical format

The fourth Top Option character `B` enables biological attenuation:

```text
'CVWB'
2
10.0 30.0 1000.0 5.0 0.04
40.0 60.0 1200.0 4.0 0.02
```

The records mean:

```text
NBioLayers
Z1 Z2 f0 Q a0
... one line per biological layer
```

Units are:

- `Z1`, `Z2`: metres;
- `f0`: hertz;
- `Q`: dimensionless;
- `a0`: dB/km.

### 5.2 Required `env_in_out.cpp` changes

`src/module/env_in_out.cpp` must:

1. Add an explicit `case 'B'` in the fourth-character switch.
2. Set `absModel` to `OceanAbsorptionModel::Biological`.
3. Consume the layer count and layer rows immediately after parsing the Top Option.
4. Complete this consumption before reading an `A` top-halfspace row or an internal-medium header.
5. Parse each biological row as exactly five numeric values.
6. Parse into a temporary `Atten_Mode`.
7. Validate the complete temporary configuration.
8. Commit it to `params.AttenUnit` only after successful parsing and validation.
9. Include the biological-layer index and offending field in error text.
10. Return `false` through the existing ENV exception path on failure.
11. Compare complete attenuation configurations when combining multiple profiles.

This file currently supplies ENV input only; no ENV writer is added by this feature.

### 5.3 ENV compatibility behavior

- A missing fourth character continues to mean `None`.
- `T` and `F` retain their existing OOK behavior.
- `B` consumes its required records.
- An unknown non-blank fourth character becomes an explicit parse failure.
- `NBioLayers=0` is valid and produces no biological contribution.
- Layer counts greater than 200 are rejected.
- Missing, incomplete, or extra-invalid biological records are rejected.

Correct handling must be verified for both:

- Biological plus a normal top boundary;
- Biological plus an `A` top halfspace.

## 6. JSON Transport

### 6.1 Canonical schema

The canonical JSON representation is:

```json
{
  "AttenUnit": {
    "AttenuationUnit": "dB/lambda",
    "OceanAbsorptionModel": "Biological",
    "BiologicalLayers": [
      {
        "Z1": 10.0,
        "Z2": 30.0,
        "f0": 1000.0,
        "Q": 5.0,
        "a0": 0.04
      }
    ]
  }
}
```

The existing `"OceanAbsorptionModel"` key is retained for backward compatibility.

### 6.2 Required `json_in_out.hpp` changes

`src/module/json_in_out.hpp` must:

1. Map `OceanAbsorptionModel::Biological` to and from `"Biological"`.
2. Add `to_json/from_json` for `BiologicalAttenuationLayer`.
3. Extend `Atten_Mode` conversion with `"BiologicalLayers"`.
4. Require the array key when the model is Biological.
5. Allow an explicitly supplied empty array.
6. Reject a non-Biological model with a non-empty layer array.
7. Derive the count from the array length.
8. Apply the same structural validation used by ENV and C++ setters.
9. Preserve the existing top-level `OOK_parameters` schema; no parallel top-level biological key is added.

### 6.3 JSON compatibility behavior

- Old non-Biological JSON remains readable.
- Non-Biological output omits `"BiologicalLayers"` so existing canonical JSON and snapshots do not change.
- New Biological JSON is fail-closed when read by an old OOK build because `"Biological"` is unknown.
- New JSON does not store a duplicate layer count.
- ENV `B` → parameters → JSON → parameters must preserve every field and layer order.

## 7. Public Interfaces

### 7.1 C++

The existing setter remains:

```cpp
void Interface::set_AttenUnit(Atten_Mode mode);
```

It validates the supplied object before changing active parameters. Failure throws `std::invalid_argument` and leaves the previous configuration unchanged.

Internal hot-path functions accept `const Atten_Mode&`.

### 7.2 Python

The pybind11 layer must expose:

- `OceanAbsorptionModel.Biological`;
- `BiologicalAttenuationLayer`;
- all five layer fields;
- `Atten_Mode.biologicalLayers`.

The `.pyi` file must describe these additions.

The high-level Python wrapper must provide a typed convenience entry point:

```python
ook_set_biological_attenuation(
    layers,
    attenuation_unit=...
)
```

Invalid input must surface as `ValueError`.

### 7.3 MATLAB

`krakenDataModel` must provide a Biological configuration helper that produces the canonical JSON structure and preserves layer order.

MATLAB-side checks should provide early feedback, while the C++ JSON import remains the authoritative validation boundary.

### 7.4 Documentation

The main README and wrapper READMEs must show:

- ENV `B` ordering;
- units for all five fields;
- JSON representation;
- Python native and high-level examples;
- MATLAB example;
- closed-interval and overlap semantics;
- the fact that `a0` is not the resonance peak because `a(f0)=a0·Q²`.

## 8. Numerical Semantics

### 8.1 Base attenuation

OOK first converts the user-specified material attenuation to Nepers per metre using the existing `N/M/m/F/W/Q/L` logic.

### 8.2 Biological attenuation

For each layer in input order:

\[
I_i(z)=
\begin{cases}
1,& Z_{1i}\le z\le Z_{2i}\\
0,& \text{otherwise}
\end{cases}
\]

\[
a_i(f,z)=
I_i(z)
\frac{a_{0i}}
{\left(1-\frac{f_{0i}^2}{f^2}\right)^2+\frac{1}{Q_i^2}}
\quad[\mathrm{dB/km}]
\]

For every matching layer, execute:

\[
\alpha_T\leftarrow\alpha_T+\frac{a_i}{8685.8896}
\quad[\mathrm{Np/m}]
\]

The operation order is normative:

1. iterate in input order;
2. perform the closed-interval test without epsilon;
3. calculate one layer's denominator;
4. calculate that layer's `a`;
5. divide that layer's value by `8685.8896`;
6. add it to the current `alphaT`;
7. continue to the next layer.

The implementation must not:

- sort, merge, or deduplicate layers;
- pre-sum all dB values before conversion;
- replace division with a precomputed reciprocal when doing so changes rounding;
- algebraically rewrite the resonance expression into an `f⁴` form;
- treat `a0` as the peak attenuation.

At resonance:

\[
a(f_0)=a_0Q^2
\]

### 8.3 Complex sound speed

Keep the current Fortran-compatible weak-loss conversion:

\[
\omega=2\pi f
\]

\[
c_I=\frac{\alpha_Tc^2}{\omega}
\]

\[
\widetilde c=c+i\,c_I
\]

No exact inverse-wavenumber reformulation is introduced.

### 8.4 SSP and halfspaces

- Biological attenuation is evaluated at stored SSP control points.
- The same configuration is applied to P-wave and S-wave calls to `CRCI`.
- Existing linear, PCHIP, and spline coefficient generation occurs after complex sound speed is updated.
- No `Z1/Z2` control points are inserted automatically.
- Overlapping layers accumulate.
- Shared endpoints are counted in both adjacent layers.
- `UpdateHSLoss()` excludes Biological attenuation.
- The existing `1e8` halfspace depth sentinel is replaced with `std::numeric_limits<double>::max()`.
- The final field continues to obtain loss from complex modal wavenumbers and `exp(-ikr)`; no extra final-pressure multiplier is added.

## 9. Validation and Failure Semantics

### 9.1 Structural rules

Every layer must satisfy:

```text
Z1, Z2, f0, Q, and a0 are finite
Z1 <= Z2
f0 > 0
Q > 0
a0 >= 0
```

Model-level rules:

- `biologicalLayers.size() <= 200`;
- Biological plus an empty list is valid;
- a non-Biological model plus a non-empty list is invalid;
- overlap is valid;
- shared endpoints are valid;
- a zero-thickness layer with `Z1==Z2` is valid;
- `a0==0` is valid;
- layer order is preserved.

### 9.2 Runtime rules

Before updating SSP loss:

- frequency must be finite and greater than zero;
- every denominator must be finite and positive;
- every layer contribution must be finite;
- cumulative attenuation must be finite;
- the resulting complex sound speed must be finite.

If the converted imaginary sound speed exceeds the real sound speed, OOK must preserve the Fortran fatal-error meaning and fail the current solve.

Fluid `cS==0` remains valid.

### 9.3 Shared validation and atomic mutation

ENV, JSON, the public setter, and runtime computation must use shared validation functions rather than duplicating rules.

Each mutating input path follows:

```text
parse into temporary object
    → validate the complete object
    → assign active state once
```

No failure may leave a partially populated Biological configuration active.

### 9.4 Error mapping

- ENV: diagnostic plus `from_env()==false`;
- JSON: diagnostic plus `from_json()==false`;
- C++ setter: `std::invalid_argument`;
- Python: `ValueError`;
- non-finite calculation: `std::domain_error` or `std::overflow_error`;
- MATLAB: early high-level diagnostic, backed by the authoritative C++ failure.

## 10. Testing Strategy

### 10.1 Formula unit tests

Use:

```text
Z1=20 m
Z2=40 m
f=f0=1000 Hz
Q=5
a0=0.04 dB/km
c=1500 m/s
base attenuation=0
```

Expected resonance attenuation:

\[
a_B=1.0\ \mathrm{dB/km}
\]

\[
\alpha_B=\frac{1}{8685.8896}
=1.1512925515\times10^{-4}\ \mathrm{Np/m}
\]

The tests must cover:

- lower boundary, interior, and upper boundary;
- points immediately outside the layer;
- overlap accumulation;
- base plus Biological addition;
- empty layer list;
- no effect for other models;
- no cross-test global state;
- SSP inclusion and halfspace exclusion.

### 10.2 ENV tests

Add positive fixtures for:

- Biological plus normal top boundary;
- Biological plus `A` top halfspace.

Add negative fixtures for:

- missing layer count;
- count/row mismatch;
- count greater than 200;
- incomplete row;
- non-numeric field;
- invalid field values;
- unknown non-blank fourth character;
- inconsistent multi-profile attenuation configuration.

### 10.3 JSON tests

Cover:

- Biological enum conversion;
- layer object conversion;
- complete `Atten_Mode` round trip;
- ENV `B` → JSON → parameters round trip;
- empty Biological array;
- missing Biological array;
- wrong field types;
- non-Biological plus non-empty array;
- over-limit array;
- unchanged output for None, Thorpe, and FrancGarr.

### 10.4 Python and MATLAB tests

Python automated tests must verify:

- native object construction;
- vector assignment;
- high-level convenience setter;
- JSON round trip;
- `ValueError` mapping.

MATLAB tests must verify:

- canonical JSON generation;
- field and array preservation;
- OOK CLI acceptance;
- invalid input rejection.

### 10.5 Fortran modal differential test

Add fixed cases:

```text
bio_uniform_off
bio_uniform_resonance
```

The Biological-on case covers the full water layer to isolate the attenuation formula from SSP boundary interpolation.

Use `kraken.exe` as the primary truth implementation.

Acceptance:

- identical mode count;
- `|ΔRe(k)| <= 1e-7 m^-1`;
- `|ΔIm(k)| <= max(1e-9, 1e-3·|Im(k_Fortran)|)`;
- Biological-on `Im(k) < 0`;
- Biological-on and Biological-off imaginary parts differ materially;
- Fortran output contains `Biological attenuation`;
- neither implementation reports fatal errors or no modes.

The new comparison tool must fail with a non-zero exit code when a threshold is exceeded.

### 10.6 End-to-end TL acceptance

Compare Biological-on versus Biological-off at 1, 5, and 10 km.

Acceptance:

- each implementation's attenuation-slope error is at most `0.02 dB/km`;
- OOK versus Fortran `P95 |ΔTL| <= 0.20 dB`;
- single-mode `max |ΔTL| <= 0.50 dB`;
- using the same Fortran `field.exe` with both MOD files gives `P95 <= 0.05 dB`;
- every pressure sample is finite and non-zero.

### 10.7 Regression gate

Completion requires:

```text
formula unit tests
+ env_in_out tests
+ json_in_out tests
+ C++/Python/MATLAB interface tests
+ Fortran modal differential test
+ TL acceptance test
+ all pre-existing CTest tests
```

## 11. Backward Compatibility and Build Impact

### 11.1 Source and data compatibility

- Existing ENV without `B` remains supported.
- Existing non-Biological JSON remains supported.
- MOD and SHD formats do not change.
- Existing C++ source using the retained fields remains source-compatible in normal use.
- Unknown non-blank ENV fourth characters change from silent fallback to explicit rejection.

### 11.2 Binary compatibility

Adding `std::vector` changes the public `Atten_Mode` object layout and the containing `OOK_parameters` layout.

This binary ABI change is accepted because there are no external C++ binary clients that must remain compatible.

The implementation release must rebuild together:

- the OOK core library;
- OOK executables;
- native tests and runners;
- the Python `.pyd`;
- any in-repository native consumers.

Pure Python files, MATLAB `.m` files, ENV files, JSON files, MOD files, and SHD files do not require compilation.

A clean build directory is recommended so stale pre-change objects cannot be mixed with new objects.

## 12. Non-Goals

This feature does not:

- modify `krakenFortran`;
- redesign all absorption models;
- correct the current Francois-Garrison environmental-parameter handling;
- change the weak-loss sound-speed approximation;
- insert biological-layer boundaries into SSP grids;
- support per-profile Biological configurations;
- add a new ENV writer;
- change MOD or SHD formats;
- implement scattering or biological sound-speed dispersion beyond the Fortran attenuation expression.

## 13. Component Map

Expected implementation areas:

- `include/OpenOceanKrakenParams.h`
- `src/algorithm/AttenMod.h`
- `src/algorithm/AttenMod.cpp`
- `src/algorithm/sspMod.h`
- `src/algorithm/sspMod.cpp`
- `src/module/env_in_out.cpp`
- `src/module/json_in_out.hpp`
- `src/module/OpenOceanKrakenInterface.cpp`
- `src/bind/OpenOceanKrakenBind.cpp`
- `wrappers/py_kraken/py_kraken/OpenOceanKraken.pyi`
- `wrappers/py_kraken/py_kraken/ook_data_model.py`
- `wrappers/py_kraken/py_kraken/ook_interface.py`
- `wrappers/py_kraken/py_kraken/__init__.py`
- `wrappers/py_kraken/tests/`
- `wrappers/m_kraken/krakenDataModel.m`
- `wrappers/m_kraken/tests/`
- `for_test/test.cpp`
- `test_for_lcov/option_examples/`
- `test_for_lcov/bad_examples/`
- `test/`
- `tools/test_bio_kraken_equivalence.py`
- `README.md`
- `wrappers/py_kraken/README.md`
- `wrappers/m_kraken/README.md`

## 14. Definition of Done

The feature is complete only when:

1. OOK can represent Biological attenuation in its public parameter model.
2. ENV `B` records are consumed in the Fortran-compatible position.
3. JSON round trips the complete model and layer list.
4. C++, Python, and MATLAB can configure the feature.
5. OOK calculates the Fortran formula with the approved operation order.
6. Halfspaces exclude Biological attenuation.
7. Invalid inputs fail atomically and consistently.
8. Multi-profile mismatches are rejected.
9. Formula, modal, and TL acceptance gates pass.
10. Existing OOK tests remain green.

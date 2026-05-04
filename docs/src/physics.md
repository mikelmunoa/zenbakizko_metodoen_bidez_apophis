# Physics

## Equations of motion

The full equations of motion for a test particle (Apophis) in an N-body system:

```math
\ddot{\mathbf{r}}_i = \sum_{j \neq i} \frac{\mu_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3}
+ \mathbf{a}_\text{PN} + \mathbf{a}_{J_2,\oplus} + \mathbf{a}_{J_3,\oplus} + \mathbf{a}_{J_4,\oplus} + \mathbf{a}_{J_2,\odot}
```

All coordinates are in the barycentric ICRF frame (J2000), in units of km and km/s.

## Post-Newtonian correction (Damour-Deruelle 1985)

The leading-order general relativistic correction from the Sun follows
Damour & Deruelle (1985):

```math
\mathbf{a}_\text{PN} = \frac{GM_\odot}{c^2 r^3}
\left[\left(\frac{4GM_\odot}{r} - v^2\right)\mathbf{r} + 4(\mathbf{r}\cdot\mathbf{v})\,\mathbf{v}\right]
```

where ``\mathbf{r}`` and ``\mathbf{v}`` are the heliocentric position and velocity,
``r = |\mathbf{r}|``, and ``c`` is the speed of light.

**Magnitude**: At 1 AU, ``|\mathbf{a}_\text{PN}| \approx 1.5 \times 10^{-13}`` km/s²,
or ~1/1220 of the Newtonian solar gravity.

**Cumulative effect on Apophis**: ~6 km position displacement by the April 2029 flyby,
~74 km/year.

## Zonal harmonics

### Earth (J₂, J₃, J₄)

The perturbation from Earth's oblateness is computed in the Earth pole frame.
The rotation uses right ascension RA = 0° and declination Dec = 90° (J2000 equator
approximation; true values: RA = 359.87°, Dec = 89.89° — difference < 0.1°).

Constants used:
- ``J_2 = 1.08263 \times 10^{-3}``
- ``J_3 = -2.5321 \times 10^{-6}``
- ``J_4 = -1.6109 \times 10^{-6}``
- ``R_\oplus = 6378.137`` km

### Sun (J₂)

- ``J_2 = 2.2 \times 10^{-7}``
- ``R_\odot = 696000`` km

The Sun's pole is at RA = 286.13°, Dec = 63.87° in ICRF.

## EMB correction

Planetary ephemerides give the Earth-Moon Barycenter (EMB, NAIF 3) position.
The true Earth center is recovered as:

```math
\mathbf{r}_\oplus = \mathbf{r}_\text{EMB} - \frac{\mu_\text{Moon}}{\mu_\text{Earth}} \mathbf{r}_\text{Moon/EMB}
```

## Physical constants

| Quantity | Value | Units |
|---|---|---|
| ``GM_\odot`` | 1.32712440018 × 10¹¹ | km³/s² |
| ``GM_\oplus`` | 398600.435 | km³/s² |
| ``c`` | 299792.458 | km/s |
| ``J_2^\oplus`` | 1.08263 × 10⁻³ | — |
| ``R_\oplus`` | 6378.137 | km |

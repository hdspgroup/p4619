# Physics-Informed Forward Soil EC Model

The current forward model predicts electrical conductivity from soil-health state variables and soil context variables.

## Response

$$
y = \log_{10}(\mathrm{EC})
$$

## Model Form

$$
\hat{y}
= \beta_0
+ \sum_i \beta_i^{(s)} \, g_i(s_i)
+ \sum_j \beta_j^{(c)} \, h_j(c_j)
+ \sum_{(i,j)\in\mathcal{I}} \beta_{ij}^{(sc)} \, g_i(s_i)\,h_j(c_j)
$$

$$
\widehat{\mathrm{EC}} = 10^{\hat{y}}
$$

## State Variables

$$
\mathbf{s} =
\left\{
N,\;P,\;K,\;\mathrm{CEC},\;\mathrm{OC},\;\mathrm{Moisture}
\right\}
$$

At present, moisture is included only if available in the data.

## Context Variables

$$
\mathbf{c} =
\left\{
\mathrm{Clay},\;\mathrm{Silt},\;\mathrm{Sand},\;\mathrm{Coarse},\;
\mathrm{pH}_{\mathrm{CaCl_2}},\;\mathrm{pH}_{\mathrm{H_2O}},\;\mathrm{CaCO_3}
\right\}
$$

## Transforms

For skewed positive variables:

$$
g(x)=\log(1+x), \qquad h(x)=\log(1+x)
$$

Otherwise the model uses the raw value:

$$
g(x)=x \qquad \text{or} \qquad h(x)=x
$$

## Selected Interaction Set

Each state variable is currently interacted with:

$$
\mathcal{I} =
\left\{
\mathrm{Clay},\;\mathrm{Silt},\;\mathrm{Sand},\;
\mathrm{pH}_{\mathrm{H_2O}},\;\mathrm{pH}_{\mathrm{CaCl_2}},\;\mathrm{CaCO_3}
\right\}
$$

Examples:

$$
g(N)\,h(\mathrm{Clay}),\quad
g(K)\,h(\mathrm{pH}_{\mathrm{H_2O}}),\quad
g(\mathrm{OC})\,h(\mathrm{CaCO_3})
$$

## Interpretation

- State variables represent soil-health quantities of interest.
- Context variables condition how those state variables translate into conductivity.
- Clay is treated as a context variable because it modifies the EC response rather than acting only as a health target.

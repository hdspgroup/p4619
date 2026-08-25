function bounds = default_property_bounds()
%DEFAULT_PROPERTY_BOUNDS Broad physical bounds for the demo inversion.
%
% These are intentionally conservative guardrails for synthetic tests.
% Keep the lower bounds physically plausible. If weakly constrained cells
% are allowed to hit zero, the property inversion can dump deeper cells to
% zero even when the synthetic prior/actual fields are nonzero everywhere.

    bounds = struct();
    bounds.N = [20, 200];
    bounds.P = [1, 150];
    bounds.K = [20, 600];
    bounds.OC = [0.1, 20];
    bounds.Moisture = [1, 60];
    bounds.Clay = [1, 100];
    bounds.Sand = [1, 100];
    bounds.Silt = [1, 100];
end

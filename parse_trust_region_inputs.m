function options = parse_trust_region_inputs(options, default_options)

options = options;

option_names = fieldnames(default_options);

for k = 1:length(option_names)
    if ~isfield(options, option_names{k})
        options.(option_names{k}) = default_options.(option_names{k});
    end
end

if options.tol_radius <= 0
    error('cmg:tol_radius', 'Option out of range: %s\n', 'tol_radius');
end
if options.eta_0 < 0 || options.eta_0 > options.eta_1
   error('cmg:eta_0', 'Option out of range: %s\n', 'eta_0'); 
end
if options.eta_1 == 0 || options.eta_1 > 1
   error('cmg:eta_1', 'Option out of range: %s\n', 'eta_1'); 
end
if options.gamma_dec <= 0 || options.gamma_dec >= 1
    error('cmg:gamma_dec', 'Option out of range: %s\n', 'gamma_dec');    
end
if options.gamma_inc <= 1
    error('cmg:gamma_inc', 'Option out of range: %s\n', 'gamma_inc');    
end
if options.eps_c <= 0
    error('cmg:eps_c', 'Option out of range: %s\n', 'eps_c');    
end
if options.criticality_mu <= options.criticality_beta
    error('cmg:criticality_mu', 'Option out of range: %s\n', ...
          'criticality_mu');    
end
if options.criticality_beta <= 0
    error('cmg:criticality_beta', 'Option out of range: %s\n', ...
          'criticality_beta');    
end
if options.criticality_omega <= 0 || options.criticality_omega >= 1
    error('cmg:criticality_omega', 'Option out of range: %s\n', ...
          'criticality_omega');    
end
if options.tol_f >= options.eps_c
    error('cmg:tol_f', 'Option out of range: %s\n', 'tol_f');    
end
 

end
function network_configs = generate_network_dataset(num_networks, varargin)
%GENERATE_NETWORK_DATASET Generate multiple random networks for dataset creation
%
% Usage:
%   configs = generate_network_dataset(100)  % 100 random networks
%   configs = generate_network_dataset(50, 'fault_probability', 0.3)
%
% Input:
%   num_networks - Number of network configurations to generate
%   varargin     - Parameters passed to generate_random_network
%
% Output:
%   network_configs - Cell array of network configuration structures

fprintf('Generating %d random network configurations...\n', num_networks);

network_configs = cell(num_networks, 1);

for i = 1:num_networks
    % Use different seed for each network
    network_configs{i} = generate_random_network(varargin{:}, 'seed', i);
    
    if mod(i, 10) == 0
        fprintf('  Generated %d/%d networks\n', i, num_networks);
    end
end

fprintf('Dataset generation complete!\n');

end 
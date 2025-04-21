clear all;
close all;
clc;

% message sequence to be encoded
message_sequence = [1 0 1 1 0 0 0];
disp('Original message sequence:');
disp(message_sequence);

constraint_length = 3;
code_generators = [6, 7];

% trellis structure
trellis_structure = poly2trellis(constraint_length, code_generators);

% Encoding
encoded_output = convenc(message_sequence, trellis_structure);
disp('Encoded bit sequence:');
disp(encoded_output);

%using encoded output as received sequence
received_sequence = encoded_output;
disp('Received sequence (no errors):');
disp(received_sequence);

% bit info
msg_length = length(message_sequence);
encoded_length = length(encoded_output);
fprintf('Message length: %d bits\n', msg_length);
fprintf('Encoded length: %d bits\n', encoded_length);

% built-in Viterbi decoder for verification
decoded_builtin = vitdec(received_sequence, trellis_structure, msg_length, 'trunc', 'hard');
disp('Decoded sequence (using built-in function):');
disp(decoded_builtin);

if isequal(message_sequence, decoded_builtin)
    disp('✓ VERIFICATION PASSED: Decoded sequence matches original message');
else
    disp('✗ VERIFICATION FAILED: Decoded sequence does not match original message');
    disp('Errors at positions:');
    disp(find(message_sequence ~= decoded_builtin));
end

%% Custom Viterbi Decoder Implementation

% decoder parameters
num_states = 2^(constraint_length-1);
time_steps = msg_length;
bits_per_symbol = 2;  % For rate 1/2 code

% received bits to symbols (pairs of bits)
rx_symbols = reshape(received_sequence, bits_per_symbol, [])';

% Viterbi algorithm data structures
path_metrics = Inf(time_steps+1, num_states);
path_metrics(1,1) = 0;  % Start at state 0 with metric 0
survivor_paths = zeros(time_steps, num_states);
branch_metrics = cell(time_steps, num_states, num_states);
valid_transitions = zeros(time_steps, num_states, num_states);

% Forward pass through trellis
disp('');
disp('FORWARD PASS - VITERBI ALGORITHM');
disp('===============================');

for t = 1:time_steps
    disp(['Time Step ' num2str(t) ':']);
    disp('Current    Input    Next     Output    Branch    Current    New      Selected');
    disp('State              State              Metric    Metric     Metric');
    disp('-----------------------------------------------------------------------');
    
    for next_state = 0:num_states-1
        for current_state = 0:num_states-1
            for input_bit = 0:1
                % Check if this input leads to next_state
                if trellis_structure.nextStates(current_state+1, input_bit+1) == next_state
                    % Mark valid transition
                    valid_transitions(t, current_state+1, next_state+1) = 1;
                    
                    % Get expected output for this transition
                    output_code = trellis_structure.outputs(current_state+1, input_bit+1);
                    expected_bits = dec2bin(output_code, bits_per_symbol) - '0';
                    
                    %  Hamming distance (branch metric)
                    branch_metric = sum(rx_symbols(t,:) ~= expected_bits);
                    branch_metrics{t, current_state+1, next_state+1} = branch_metric;
                    
                    %  new path metric
                    new_metric = path_metrics(t, current_state+1) + branch_metric;
                    
                    % Format output for display
                    output_str = [num2str(expected_bits(1)) num2str(expected_bits(2))];
                    
                    % Display computation
                    if ~isinf(path_metrics(t, current_state+1))
                        fprintf('  %d        %d        %d       %s        %d        %.1f       %.1f', ...
                            current_state, input_bit, next_state, output_str, branch_metric, ...
                            path_metrics(t, current_state+1), new_metric);
                        
                        % Update path metric if better
                        if new_metric < path_metrics(t+1, next_state+1)
                            path_metrics(t+1, next_state+1) = new_metric;
                            survivor_paths(t, next_state+1) = current_state;
                            fprintf('      ✓\n');
                        else
                            fprintf('\n');
                        end
                    end
                end
            end
        end
    end
    
    % path metrics after this time step
    disp(' ');
    disp(['Path Metrics after Time Step ' num2str(t) ':']);
    for s = 0:num_states-1
        if ~isinf(path_metrics(t+1, s+1))
            fprintf('  State %d: %.1f\n', s, path_metrics(t+1, s+1));
        else
            fprintf('  State %d: Inf\n', s);
        end
    end
    disp(' ');
end

% Traceback to find survivor path
disp('TRACEBACK - FINDING SURVIVOR PATH');
disp('================================');

%  state with minimum path metric at final time step
[~, best_state_idx] = min(path_metrics(time_steps+1, :));
final_state = best_state_idx - 1;

fprintf('Final state with minimum metric: State %d (Metric = %.1f)\n', ...
    final_state, path_metrics(time_steps+1, best_state_idx));

% Initializing traceback variables
state_sequence = zeros(1, time_steps+1);
state_sequence(end) = final_state;
decoded_output = zeros(1, time_steps);

% Performing traceback
disp(' ');
disp('Traceback Process:');
disp('Time    Current    Previous    Input');
disp('Step    State      State       Bit');
disp('--------------------------------');

for t = time_steps:-1:1
    current_state = state_sequence(t+1);
    previous_state = survivor_paths(t, current_state+1);
    state_sequence(t) = previous_state;
    
    % input bit that caused this transition
    for input_bit = 0:1
        if trellis_structure.nextStates(previous_state+1, input_bit+1) == current_state
            decoded_output(t) = input_bit;
            break;
        end
    end
    
    fprintf('  %d      %d          %d          %d\n', ...
        t, current_state, previous_state, decoded_output(t));
end

% Displaying custom decoder results
disp(' ');
disp('Custom decoder output:');
disp(decoded_output);

% Verifying custom decoder
if isequal(message_sequence, decoded_output)
    disp('✓ CUSTOM DECODER VERIFICATION PASSED');
else
    disp('✗ CUSTOM DECODER VERIFICATION FAILED');
    disp('Errors at positions:');
    disp(find(message_sequence ~= decoded_output));
end

%% Visualization of Trellis Diagram

% figure for trellis diagram
figure('Name', 'Trellis Diagram', 'Position', [100, 100, 900, 600]);
hold on;
grid off;
box on;

% Spacing parameters
time_spacing = 1;
state_spacing = 1.2;
node_size = 8;

% Plot all states at each time step
for t = 1:time_steps+1
    for s = 0:num_states-1
        % Plot node
        plot(t*time_spacing, -s*state_spacing, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', node_size);
        
        % Add path metric label
        if ~isinf(path_metrics(t, s+1))
            text(t*time_spacing, -s*state_spacing-0.25, sprintf('PM: %.1f', path_metrics(t, s+1)), ...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', 'blue');
        end
    end
end

% Plot all transitions
for t = 1:time_steps
    for from_state = 0:num_states-1
        for to_state = 0:num_states-1
            if valid_transitions(t, from_state+1, to_state+1)
                % Find input bit for this transition
                for input_bit = 0:1
                    if trellis_structure.nextStates(from_state+1, input_bit+1) == to_state
                        % Get output bits
                        output_code = trellis_structure.outputs(from_state+1, input_bit+1);
                        output_bits = dec2bin(output_code, bits_per_symbol);
                        
                        % Get branch metric
                        bm = branch_metrics{t, from_state+1, to_state+1};
                        
                        % Plot regular transition (gray dashed)
                        plot([t, t+1]*time_spacing, [-from_state, -to_state]*state_spacing, ...
                            'k--', 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7]);
                        
                        % Add transition label
                        mid_x = (t + 0.5) * time_spacing;
                        mid_y = (-from_state - 0.5 * (to_state - from_state)) * state_spacing;
                        
                        label = sprintf('%d/%s (BM:%d)', input_bit, output_bits, bm);
                        text(mid_x, mid_y+0.1, label, 'FontSize', 7, 'Color', 'black', ...
                            'HorizontalAlignment', 'center');
                        
                        break;
                    end
                end
            end
        end
    end
end

% Highlight survivor path
for t = 1:time_steps
    current_state = state_sequence(t);
    next_state = state_sequence(t+1);
    
    % Plot survivor path with thick green line
    plot([t, t+1]*time_spacing, [-current_state, -next_state]*state_spacing, ...
        'g-', 'LineWidth', 2.5);
    
    % Highlight nodes on the path
    plot(t*time_spacing, -current_state*state_spacing, 'go', ...
        'MarkerFaceColor', 'g', 'MarkerSize', node_size+2);
end
% Highlight final node
plot((time_steps+1)*time_spacing, -state_sequence(end)*state_spacing, 'go', ...
    'MarkerFaceColor', 'g', 'MarkerSize', node_size+2);

% Add annotations and formatting
title('Trellis Diagram with Path Metrics and Survivor Path');
xlabel('Time Step');
ylabel('State');
set(gca, 'YTick', -state_spacing * (0:num_states-1));
set(gca, 'YTickLabel', 0:num_states-1);
set(gca, 'XTick', 1:time_steps+1);
xlim([0.5, time_steps+1.5]);
ylim([-(num_states-0.5)*state_spacing, 0.5]);

%% Test with specific received sequence from assignment

% Given received sequence from assignment
specific_received = [0 1 1 1 0 1 1 1 0 1 0 1 1 1];
disp(' ');
disp('TESTING WITH SPECIFIC RECEIVED SEQUENCE');
disp('======================================');
disp('Assignment received bits:');
disp(specific_received);

% Validate length
expected_length = 2 * msg_length;
if length(specific_received) ~= expected_length
    fprintf('WARNING: Expected %d bits for %d input bits, but got %d bits.\n', ...
        expected_length, msg_length, length(specific_received));
end

% Decode with built-in decoder
decoded_specific = vitdec(specific_received, trellis_structure, msg_length, 'trunc', 'hard');
disp('Decoded result from specific sequence:');
disp(decoded_specific);

% Verify against original message
if isequal(message_sequence, decoded_specific)
    disp('✓ Specific sequence decodes to original message');
else
    disp('✗ Specific sequence does not decode to original message');
    disp('Differences at positions:');
    disp(find(message_sequence ~= decoded_specific));
end

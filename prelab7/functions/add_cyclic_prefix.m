function ofdm_symbol_cp = add_cyclic_prefix(ofdm_symbol, cp_length)
    % Add a cyclic prefix to an OFDM symbol block
    
    % Inputs:
    % - ofdm_symbol: The OFDM symbol block without cyclic prefix
    % - cp_length: Length of the cyclic prefix (in samples)
    
    % Output:
    % - ofdm_symbol_cp: The OFDM symbol block with cyclic prefix
    
    % Get the length of the OFDM symbol block
    
    % Add cyclic prefix
    cyclic_prefix = ofdm_symbol(end - cp_length + 1:end); % Copy the last 'cp_length' samples
    ofdm_symbol_cp = [cyclic_prefix; ofdm_symbol];        % Append cyclic prefix to the beginning
end

function ofdm_symbol_block = remove_cyclic_prefix(ofdm_symbol, cp_length)
    % Remove cyclic prefix
    ofdm_symbol_block = ofdm_symbol((cp_length+1):end); % Copy the last 'cp_length' samples
end

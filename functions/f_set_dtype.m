function Y = f_set_dtype(Y, dtype)

if strcmpi(dtype, 'uint8')
    Y = uint8(Y);
elseif strcmpi(dtype, 'uint16')
    Y = uint16(Y);
elseif strcmpi(dtype, 'int8')
    Y = int8(Y);
elseif strcmpi(dtype, 'int16')
    Y = int16(Y);
end

end
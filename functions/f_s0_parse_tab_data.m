function AC_data2 = f_s0_parse_tab_data(params)

AC_data = readtable(params.dset_table_fpath);

%%
idx_idx = strcmpi(AC_data.Properties.VariableNames, 'Idx');
if sum(idx_idx)
    idx1 = ~isnan(AC_data.(AC_data.Properties.VariableNames{idx_idx}));
    AC_data = AC_data(idx1,:);
end

idx_proc = strcmpi(AC_data.Properties.VariableNames, 'do_proc');
if sum(idx_proc)
    idx1 = AC_data.(AC_data.Properties.VariableNames{idx_proc}) == 1;
    AC_data = AC_data(idx1,:);
end

%%
fields1 = fields(params.limit);

AC_data2 = AC_data;
for n_fl = 1:numel(fields1)
    field_cur = fields1{n_fl};
    tag1 = params.limit.(field_cur);
    idx_fl = strcmpi(AC_data.Properties.VariableNames, field_cur);
    if sum(idx_fl)
        field_cur2 = AC_data.Properties.VariableNames{idx_fl};
        if numel(tag1)
            if isnumeric(tag1)
                idx1 = AC_data2.(field_cur2) == tag1;
                
            else
                idx1 = strcmpi(AC_data2.(field_cur2), tag1);
            end
            AC_data2 = AC_data2(idx1,:);
        end
    end
end


end
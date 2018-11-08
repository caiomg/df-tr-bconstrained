function result = is_lambda_poised(model, options)
    % IS_LAMBDA_POISED tests wether a model is lambda-poised for the
    % given options
    %
    % Important: it assumes there is a finite radius_max defined.

    % pivot_threshold defines how well-poised we demmand a model to be
    pivot_threshold = options.pivot_threshold;
    [dim, points_num] = size(model.points_abs);

    if strcmp(options.basis, 'dummy')
        % If not modeling...
        result = true;
    else
        if points_num >= dim + 1
            % Fully linear, already
            result = true;
%             % But let's double-check
%             if ~isempty(find(model.pivot_absvalues(1:points_num) < pivot_threshold, 1))
%                 warning('cmg:pivot_not_satisfied', 'Low pivot');
%             end
        else
            result = false;
        end
    end
end

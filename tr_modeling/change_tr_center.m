function [model, exitflag] = change_tr_center(model, new_point, ...
                                              new_fvalues, options)
% CHANGE_TR_CENTER - accepts a step of the algorithm
% Calls the appropriate functions for model maintenance so the
% trust-region is shifted to this new center

    STATUS_POINT_ADDED = 1;
    STATUS_POINT_EXCHANGED = 2;
    STATUS_MODEL_REBUILT = 4;

    model_is_complete = is_complete(model);
    if ~model_is_complete
        % Add this point
        relative_pivot_threshold = options.pivot_threshold;
        [model, point_added] = add_point(model, new_point, new_fvalues, ...
                                         relative_pivot_threshold);
        if point_added
            % Function add_point adds this as the last.
            % Now we have to set as TR center
            model.tr_center = size(model.points_abs, 2); % Last among
                                                         % the points
            exitflag = STATUS_POINT_ADDED;
        end
    end
    if model_is_complete || ~point_added
        % Try to exchange a point of the model with this new one
        relative_pivot_threshold = options.pivot_threshold;
        [model, point_exchanged, pt_i] = exchange_point(model, ...
                                                        new_point, ...
                                                        new_fvalues, ...
                                                        relative_pivot_threshold,...
                                                        true);
        if point_exchanged
            % Make new point the center
            model.tr_center = pt_i;
            exitflag = STATUS_POINT_EXCHANGED;
        else
            % add_point and exchange_point failed
            % We still need to add this new point as TR center
            % Model needs rebuilding
            model.points_abs(:, end+1) = new_point;
            model.fvalues(:, end+1) = new_fvalues;
            model.tr_center = size(model.points_abs, 2); % Last
            model = rebuild_model(model, options);
            exitflag = STATUS_MODEL_REBUILT;
        end
    end

end
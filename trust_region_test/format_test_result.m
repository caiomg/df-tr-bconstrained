function s = format_test_result(test_result)

    s = sprintf('| %10s |', test_result.name);
    if isfield(test_result.test, 'fx')
        s = [s, sprintf(' Rf(x): % +9.3g |', test_result.ref.fx)];
        s = [s, sprintf(' tf(x): % +9.3g |', test_result.test.fx)];
        
        s = [s, sprintf(' Revals: % 6d |', test_result.ref.count)];
        s = [s, sprintf(' tevals: % 6d |', test_result.test.count)];

        s = [s, sprintf(' Rviol.: % 9.3g |', test_result.ref.viol)];
        s = [s, sprintf(' tviol.: % 9.3g |', test_result.test.viol)];
    else
        s = [s, sprintf(' Rf(x):           |')];
        s = [s, sprintf(' tf(x):           |')];
        s = [s, sprintf(' Revals:        |')];
        s = [s, sprintf(' tevals:        |')];
        s = [s, sprintf(' Rviol.:           |')];
        s = [s, sprintf(' tviol.:           |')];
    end
    s = [s, newline];
end
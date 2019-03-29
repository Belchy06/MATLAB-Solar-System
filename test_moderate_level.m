function test_moderate_level(unit_under_test)
% TEST_MODERATE_LEVEL Test the simulator against the moderate level of achievement.
%
% TEST_MODERATE_LEVEL(@unit) tests a function called "unit" instead of the
% default, "solarsystem".
%
% This is provided as a means for your to test your program's accuracy. We
% supply here high precision answers that can use as a benchmark against
% which to compare your code.
%
% A program similar to this one will be used during marking to test your
% program's accuracy and speed. We'll use different initial conditions, so
% don't try simply hard-coding these answers! :)
%

% Default to a function named "solarsystem"
if nargin < 1
    unit_under_test = @solarsystem;
end

% Physical constants

% Sun
% Data from http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
data.sun.p = [0 0]; % put the sun at the origin
data.sun.v = [0 0]; % no velocity
data.sun.mass = 1988500e24; % kg

% Earth
% Data from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
%
% Explanation: when the Earth is closest to the sun ("perihelion") its
% orbital velocity is maximised. Start the simulation at perihelion
% using the data from the above webpage.
data.earth.p = [147.09e9 0];
data.earth.v = [0 30.29*1000];
data.earth.mass = 5.9723e24;

% A comet
% Data completely made up
phi = 45 * pi / 180;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)]; % 2D rotation matrix
data.comet.p = (R * data.earth.p')';
data.comet.v = 0.4 * (R * data.earth.v')';
data.comet.mass = 5e21;

% Test 1
Test_2D_Sun_Earth_Comet(false);

% Test 2
Test_2D_Sun_Earth_Comet(true);

    function test_result(parameter, value, units, comparator, benchmark)
        if strcmp(units, '%')
            fprintf('  %28s :  %-15.6f', [parameter ' (' units ')'], value);
        else
            fprintf('  %28s :  %-15.6g', [parameter ' (' units ')'], value);
        end

        if nargin == 5
            if comparator(value, benchmark)
                fprintf('   ** PASS. Meets or exceeds the expectation of %g%s', benchmark, units);
            else
                fprintf('   ** FAIL. Does not meet the expectation of %g%s', benchmark, units);
            end
        end
        fprintf('\n');
    end


    function Test_2D_Sun_Earth_Comet(speed_test)
        if speed_test
            fprintf('<strong>*** [Moderate level] Three body simulation with Sun, Earth and a comet (execution speed test)</strong>\n');
        else
            fprintf('<strong>*** [Moderate level] Three body simulation with Sun, Earth and a comet</strong>\n');
        end
        
        p = [data.sun.p; data.earth.p; data.comet.p];
        v = [data.sun.v; data.earth.v; data.comet.v];
        mass = [data.sun.mass; data.earth.mass; data.comet.mass];

        % Run the program
        tic();
        [final_p, final_v] = unit_under_test(p, v, mass, 600*24*60*60, speed_test);
        t = toc();
        test_result('Execution time', t, 's');
        
        % Check the dimensions of the return values
        assert(all(size(final_p) == [3 2]), 'Expected size of return value "p" to be 3x2, received %ix%i instead. Did you implement 3 or more celestial bodies?', size(final_p, 1), size(final_p, 2));
        assert(all(size(final_v) == [3 2]), 'Expected size of return value "v" to be 3x2, received %ix%i instead. Did you implement 3 or more celestial bodies?', size(final_p, 1), size(final_p, 2));
        
        % Check the answers
        correct_p = [735284.879935314 5063046.74917656;-98070212939.0299 -115160878452.147;71150402679.2697 110734051521.905];
        correct_v = [-0.0680778317493048 0.147475054352529;22675.1355951941 -18802.003103355;-18556.0820979253 -3827.78340659555];
        
        expectations = [0.5 0.5 30]; % the comet is much harder to simulate due to its fast motion near the sun
        for i = 1:size(p,1)
            test_result(sprintf('Object %i position error', i), norm(final_p(i,:) - correct_p(i,:))/norm(correct_p(i,:))*100, '%', @le, expectations(i));
            test_result(sprintf('Object %i velocity error', i), norm(final_v(i,:) - correct_v(i,:))/norm(correct_v(i,:))*100, '%', @le, expectations(i));
        end  
    end

end
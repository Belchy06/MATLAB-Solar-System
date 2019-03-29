function test_base_level(unit_under_test)
% TEST_BASE_LEVEL Test the simulator against the base level of achievement.
%
% TEST_BASE_LEVEL(@unit) tests a function called "unit" instead of the
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


% Test 1
fprintf('<strong>*** [Basic level] Two body simulation with Sun and Earth</strong>\n');

p = [data.sun.p; data.earth.p];
v = [data.sun.v; data.earth.v];
mass = [data.sun.mass; data.earth.mass];

% Run the program for half an orbit
tic();
[final_p, final_v] = unit_under_test(p, v, mass, 365.242*24*60*60 / 2);
t = toc();
test_result('Execution time', t, 's');

% Check the dimensions of the return values
assert(all(size(final_p) == [2 2]), 'Expected size of return value "p" to be 2x2, received %ix%i instead', size(final_p, 1), size(final_p, 2));
assert(all(size(final_v) == [2 2]), 'Expected size of return value "v" to be 2x2, received %ix%i instead', size(final_p, 1), size(final_p, 2));

% The correct answer (as computed using more sophisticated
% mathematics).
correct_p = [ 8.988499926438609e+05   1.434697390203846e+06; ...
             -1.521855237299393e+11   2.414309905830630e+08];
correct_v = [ 1.410614696972120e-04   0.178899877935263; ...
             -46.966952847804066     -2.927539478496893e+04];
% ^^ hey look, it's the answers!
% Don't try to be too clever here. Obviously we'll use different
% inputs for actual marking ;)

test_result('Positional error', norm(final_p(2,:) - correct_p(2,:))/norm(correct_p(2,:))*100, '%', @le, 0.5);
test_result('Velocity error', norm(final_v(2,:) - correct_v(2,:))/norm(correct_v(2,:))*100, '%', @le, 0.5);


% Test 2
fprintf('<strong>*** [Basic level] Two body simulation of binary stars</strong>\n');

p = [-170e9 0; 98e9 0];
v = [0 -10e3];
mass = [1.1*data.sun.mass; 0.8*data.sun.mass];
% compute the velocity of object #2 such that the total momentum is
% zero (and the simulation stays near the origin)
%        m_1 v_1 + m_2 v_2 = 0
%       -m_1 v_1 / m_2 = v_2
v(2,:) = -mass(1) / mass(2) * v(1,:);

% Run the program
tic();
[final_p, final_v] = unit_under_test(p, v, mass, 1000*24*60*60);
t = toc();
test_result('Execution time', t, 's');

% Check the dimensions of the return values
assert(all(size(final_p) == [2 2]), 'Expected size of return value "p" to be 2x2, received %ix%i instead', size(final_p, 1), size(final_p, 2));
assert(all(size(final_v) == [2 2]), 'Expected size of return value "v" to be 2x2, received %ix%i instead', size(final_p, 1), size(final_p, 2));

% The correct answer (as computed using more sophisticated
% mathematics).
correct_p = [-4.723705790100282e+10   6.291039247339576e+10; ...
             -7.079904538612178e+10  -8.650178965092014e+10];
correct_v = [-1.647419489652207e+04   9.275728200331376e+03; ...
              2.265201798271780e+04  -1.275412627545568e+04];        

test_result('Object 1 position error', norm(final_p(1,:) - correct_p(1,:))/norm(correct_p(1,:))*100, '%', @le, 10);
test_result('Object 1 velocity error', norm(final_v(1,:) - correct_v(1,:))/norm(correct_v(1,:))*100, '%', @le, 10);
test_result('Object 2 position error', norm(final_p(2,:) - correct_p(2,:))/norm(correct_p(2,:))*100, '%', @le, 10);
test_result('Object 2 velocity error', norm(final_v(2,:) - correct_v(2,:))/norm(correct_v(2,:))*100, '%', @le, 10);     


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
end

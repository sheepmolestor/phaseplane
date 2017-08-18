function phaseplane2(func, varargin)
%PHASEPLANE2   Visualise 2D (time-dependent) vector field and particle
%trajectories
%    
%    PHASEPLANE2(F)   F is a function handle to a function of two
%    variables, F(t, x), which takes in a 2xN array of points and returns a
%    2xN array of velocities
%
%Example: Make sure you've saved pendulum.m and doublegyre.m in the same
%directory, then
%    phaseplane2(@pendulum)
%    phaseplane2(@doublegyre, 'xlim', [0 2], 'ylim', [0 1], 'tlim', [0 10])

% Parse input parameters
p = inputParser();
p.addParameter('xlim', [-10 10]);
p.addParameter('ylim', [-10 10]);
p.addParameter('tlim', [0 300]);
p.addParameter('ngrid', 20);
p.addParameter('xpt0', []);
p.addParameter('speedup', 1);
p.parse(varargin{:});

if p.Results.speedup ~= 1
    fcn = @(t, x) p.Results.speedup * func(t, x);
else
    fcn = func;
end

% Define the grid
[xgrid, ygrid] = meshgrid(...
    linspace(p.Results.xlim(1), p.Results.xlim(2), p.Results.ngrid), ...
    linspace(p.Results.ylim(1), p.Results.ylim(2), p.Results.ngrid));
xgrid = xgrid(:)';
ygrid = ygrid(:)';
if ~isempty(p.Results.xpt0)
    xpt = p.Results.xpt0;
else
    xpt = zeros(2, 0);
end
% Get initial velocity field
t = p.Results.tlim(1);
V = fcn(t, [xgrid; ygrid]);

% Display the vector field
hF = figure();
hA = axes();
hQ = quiver(hA, xgrid, ygrid, V(1, :), V(2, :), 'ShowArrowHead', 'off');
hold on
hX = plot(0, 0, 'r.',  'markersize', 30);
hold off
set(hX, 'Xdata', xpt(1, :), 'Ydata', xpt(2, :));
%axis equal
set(gca, 'xlim', 1.05 * p.Results.xlim, 'ylim', 1.05 * p.Results.ylim)
grid on
set(hF, 'renderer', 'opengl');
tic;

% Set up a system timer
set(hF, 'WindowButtonDownFcn', @newpoint, 'CloseRequestFcn', @closereq)

hTimer = timer( ...
    'ExecutionMode', 'fixedspacing', ...
    'Period',        1e-3, ...
    'TimerFcn',      @next, ...
    'StartDelay',    0);
start(hTimer)


% NESTED FUNCTIONS (CALLBACKS)
    function newpoint(varargin)
        % Get current point
        cp = get(hA, 'CurrentPoint');
        xpt = [xpt cp(1, 1:2)'];
    end

    function next(varargin)
        t = toc();
        
        if t > p.Results.tlim(2)
            stop(hTimer)
            delete(hTimer)
            return
        end
        
        dt = get(hTimer, 'InstantPeriod');
        if isnan(dt)
            dt = 0.05;
        end
        
        % Update all particle positions (x is a 2xN array of particles)
        % using a fourth order RK solver
        k1 = fcn(t, xpt) * dt;
        k2 = fcn(t, xpt + 0.5*k1) * dt;
        k3 = fcn(t, xpt + 0.5*k2) * dt;
        k4 = fcn(t, xpt + k3) * dt;
        xpt = xpt + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        
        % Delete any particles that have disappeared off the grid
        idx = xpt(1, :) < p.Results.xlim(1) | ...
            xpt(1, :) > p.Results.xlim(2) | ...
            xpt(2, :) < p.Results.ylim(1) | ...
            xpt(2, :) > p.Results.ylim(2);
        xpt = xpt(:, ~idx);
        
        % Update the flow vectors
        V = fcn(t, [xgrid; ygrid]);
        
        % Update the graph
        set(hX, 'XData', xpt(1, :), 'YData', xpt(2, :));
        set(hQ, 'UData', V(1, :), 'VData', V(2, :));
        drawnow
    end

    function closereq(varargin)
        out = timerfindall;
        if ~isempty(out)
            stop(out)
            delete(out)
        end
        delete(varargin{1})
    end
end %phaseplane2
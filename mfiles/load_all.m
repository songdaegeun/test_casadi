function load_all(varargin)
    global XX YY ZZ x_opt y_opt circle_x circle_y line_x line_y;

    data1 = load("mat/"+varargin{1});  % Loads XX, YY, ZZ
    XX = data1.XX;
    YY = data1.YY;
    ZZ = data1.ZZ;
    data2 = load("mat/"+varargin{2});  % Loads x_opt, y_opt
    x_opt = data2.x_opt;
    y_opt = data2.y_opt;
    data3 = load("mat/"+varargin{3});  % Loads circle_x, circle_y
    circle_x = data3.circle_x;
    circle_y = data3.circle_y;
    data4 = load("mat/"+varargin{4});  % Loads line_x, line_y
    line_x = data4.line_x;
    line_y = data4.line_y;
end
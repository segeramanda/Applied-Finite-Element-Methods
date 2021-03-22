function f=my_f(x)
    if abs(x)<=0.1
        f = -5;

    elseif abs(0.2-abs(x))<=0.1
        f = 25;

    elseif abs(0.6-abs(x))<=0.1
        f = -30;

    elseif abs(0.9-abs(x))<=0.2
        f = 20;
    else 
        f=1;
    end
end
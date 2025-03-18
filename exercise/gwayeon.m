function gwayeon(this, that)
    disp("This " + this);
    disp("That " + that);

    if (abs(this - that) < 1E-6)
        disp("Hooray");
    else
        disp("NOO");
    end

end

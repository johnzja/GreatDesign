%% Construct S4 symmetry group using a, b, c Three elements.
S4 = perms([0,1,2,3]).';

a = [0,2,1,3].';
b = [1,3,0,2].';
c = [0,2,3,1].';

% Label each element of S4 with the form a^(n1)*b^(n2)*c^(n3).
global monomial_table_abc;
global monomial_table_cab;
monomial_table_abc = zeros(2,4,3);

% Label each element of S4 with the form c^(n1)*a^(n2)*b^(n3).
monomial_table_cab = zeros(3,2,4);

% mult_table.
mult_table = zeros(24,24);
for idx_1 = 1:24
    for idx_2 = 1:24
        g1g2 = Gmult(S4(:,idx_1), S4(:,idx_2));
        for idx = 1:24
            if all(g1g2 == S4(:,idx))
                mult_table(idx_1, idx_2) = idx;
                break;
            end
        end
    end
end

global n_power_abc;
global n_power_cab;
n_power_abc = zeros(3,24);
n_power_cab = zeros(3,24);

for idx_a = 0:1
    for idx_b = 0:3
        for idx_c = 0:2
            % Calculate the group element.
            g1 = Gmult(Gmult(Gpower(a,idx_a),Gpower(b,idx_b)),Gpower(c,idx_c));
            for idx = 1:24
                if all(g1==S4(:,idx))
                    monomial_table_abc(idx_a+1,idx_b+1,idx_c+1) = idx;
                    n_power_abc(:,idx) = [idx_a, idx_b, idx_c].';
                    break;
                end
            end
            g2 = Gmult(Gmult(Gpower(c,idx_c),Gpower(a,idx_a)),Gpower(b,idx_b));
            for idx = 1:24
                if all(g2==S4(:,idx))
                    monomial_table_cab(idx_c+1,idx_a+1,idx_b+1) = idx;
                    n_power_cab(:,idx) = [idx_c, idx_a, idx_b].';
                    break;
                end
            end
        end
    end
end


% Given the element idx, output the powers of the element.
% bc = c^2*b^3: Verified.

%% Calculate the blackwell-measure W+ channel with this S4 group.
% Form: 'abc'
% [pi_0, pi_1, pi_2, pi_3].
pi_y_prime = [0, 1, 1, 0;...
              0, 1, 3, 2;...
              0, 0, 0, 0];
          
tau_u_prime = [0, 1, 1, 0;...
               0, 1, 3, 2;...
               0, 0, 0, 0];

c_square_idx = monomial_table_abc(0+1, 0+1, 1+2);

for u1 = 0:3
    for y1 = 0:3
        for y2 = 0:3
            pi_y1 = monomial_table_abc(pi_y_prime(1,y1+1)+1, pi_y_prime(2,y1+1)+1,pi_y_prime(3,y1+1)+1);
            tau_u1 = monomial_table_abc(tau_u_prime(1, u1+1)+1,tau_u_prime(2, u1+1)+1,tau_u_prime(3, u1+1)+1);
            p_ops = mult_table(mult_table(pi_y1, tau_u1), c_square_idx);
            q_ops = monomial_table_abc(pi_y_prime(1,y2+1)+1, pi_y_prime(2,y2+1)+1,pi_y_prime(3,y2+1)+1);
            fprintf('(');
            display_element(p_ops, 'cab');
            fprintf(', ');
            display_element(q_ops, 'cab');
            fprintf(')\t');
        end
    end
    fprintf('\n');
end

%% 
function ret=Gmult(g1,g2)
    ret=zeros(4,1);
    for idx = 1:4
        ret(idx) = g2(1+g1(idx));
    end
end

function ret=Gpower(g,n)
    ret = [0,1,2,3].';
    if n == 0
        return;
    else
        
        for idx = 1:n
            ret = Gmult(ret, g);
        end
    end
end


function display_element(g_idx, method)
global n_power_abc;
global n_power_cab;
    if strcmp(method, 'abc')
        powers = n_power_abc(:,g_idx);
        if ~any(powers)
            fprintf('1');
            return;
        end
        if powers(1) > 0
            fprintf('a');
            if powers(1) > 1
                fprintf('^%d', powers(1));
            end
        end
        if powers(2) > 0
            fprintf('b');
            if powers(2) > 1
                fprintf('^%d', powers(2));
            end
        end
        if powers(3) > 0
            fprintf('c');
            if powers(3) > 1
                fprintf('^%d', powers(3));
            end
        end
    elseif strcmp(method, 'cab')
        powers = n_power_cab(:,g_idx);
        if ~any(powers)
            fprintf('1');
            return;
        end
        if powers(1) > 0
            fprintf('c');
            if powers(1) > 1
                fprintf('^%d', powers(1));
            end
        end
        if powers(2) > 0
            fprintf('a');
            if powers(2) > 1
                fprintf('^%d', powers(2));
            end
        end
        if powers(3) > 0
            fprintf('b');
            if powers(3) > 1
                fprintf('^%d', powers(3));
            end
        end
    end
end


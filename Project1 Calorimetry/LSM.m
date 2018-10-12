function [ m b sig_y sig_b sig_m ] = LSM(t,y)

b_coef = ones(1,length(t))';
A = [ t b_coef ];
b = y;


x = inv(transpose(A) * (A)) * transpose(A) * b;
m = x(1);
b = x(2);

%% sigma y is equation  from the book [ ? (1/N-2) sum(yi-mti-b) ]


% the under square root
for i=1:length(b)
    
    sq(i) = (b(i) - x(1)*A(i,1) - x(2))^2;
    
end

% sig y 

sig_y = sqrt( (1/(length(b)-2)) * sum(sq));

% weight matrix 
weight_matrix = zeros(length(b),length(b));

%put the sig/y in the daiagonal 

for i = 1:length(b)
    
    weight_matrix(i,i) = 1 / (sig_y)^2;
    
end

Q = inv( transpose(A) * weight_matrix * A ) ;

sig_m = sqrt(Q(1,1));
sig_b = sqrt(Q(2,2));


end
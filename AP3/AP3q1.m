#intervalo = [inicio, fim]
#S = [v; x]

function resultado = AP3q1(S0, intervalo, delta_t, erro, w, m)
  Sant = zeros(2, 1);
  Snovo = S0;
  ant = 0;
  novo = 1;
  k = 0;
  
  while (abs((novo - ant) / novo) >= erro && k < 10000)
    ant = novo;
    Sant = Snovo;
    
    Snovo = S0;
    for Ti = intervalo(1,1) : delta_t : intervalo(1,2)
      Snovo = range4kutta(Snovo, Ti, delta_t, w, m);
    endfor
    
    novo = Snovo(2, 1);
    delta_t /= 2;
    k++;
  endwhile
  
  resultado = Snovo;
endfunction
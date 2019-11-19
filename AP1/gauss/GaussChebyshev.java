//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
import java.lang.Math;
public class GaussChebyshev {

	public static void main(String[] args) {
		GaussChebyshev gc = new GaussChebyshev();					//GRAU
		System.out.println("\nO valor da integral é: " + gc.chebyshev(10) + "\n"); 
	}
	
	private double f(double x)
	{
	    //Usando para os testes a função x².
	    //A integral de sqrt(1-x²) * x² no
	    //intervalo [-1;1] é igual à 1.5708
	    //Note que, como essa função é de grau
	    //dois, o número de pontos necessários
	    //pra calcular essa integral é pelo
	    //menos 2.
	    return Math.sin(x);
	}


	public double chebyshev(int degree)
	{
	    double sum = 0;
	    double xk;

	    for (int i=0; i<degree; i++)
	    {
	        xk = Math.cos((2*i-1)*Math.PI/(2*degree));
	        sum = sum + f(xk);
	    }

	    sum = sum*Math.PI/degree;

	    return sum;
	}


}

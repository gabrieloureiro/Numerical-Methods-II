//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
public class NewtonCotesFechada {

	public static void main(String[] args) {
		NewtonCotesFechada ncf = new NewtonCotesFechada();						//lim inf, lim sup, erro, NCF_grau
		System.out.println("\nO valor da integral é: " + ncf.integralPartitions(1, 10, 0.001, 3) + "\n"); 
	}
	
	private double f(double x)
	{
		// Usando para os testes a função 30x⁴ - 2x³ + 8x - 12.5.
		// A integral dessa função no intervalo de -1 a 1 é -13.
		//return (30 * Math.pow(x, 4) - 2 * Math.pow(x, 3) + 8 * x - 12.5);
		return 1 / (9 * Math.pow(x, 2));
	}
	
	public double integralGrau2(double inferior, double superior)
	{
		double f1;
	    double f2;
	    double f3;
	    double value;
	    double deltaX;

	    deltaX = (superior - inferior);
	    f1 = f(inferior);
	    f2 = f((inferior + superior)/2);
	    f3 = f(superior);
	    value = (deltaX/6)*(f1 + 4*f2 + f3);

	    return value;
	}
	
	
	public double integralGrau3(double inferior, double superior)
	{
	    double f1;
	    double f2;
	    double f3;
	    double f4;
	    double value;
	    double deltaX;

	    deltaX = (superior - inferior);
	    f1 = this.f(inferior);
	    f2 = this.f(inferior + deltaX/3);
	    f3 = this.f(inferior + deltaX*2/3);
	    f4 = this.f(superior);
	    value = (deltaX*3/24)*(f1 + 3*f2 + 3*f3 + f4);

	    return value;
	}
	
	public Double integralGrau4(double inferior, double superior)
	{
	    double f1;
	    double f2;
	    double f3;
	    double f4;
	    double f5;
	    double value;
	    double deltaX;

	    deltaX = (superior - inferior);
	    f1 = f(inferior);
	    f2 = f(inferior + deltaX/4);
	    f3 = f(inferior + deltaX/2);
	    f4 = f(inferior + deltaX*3/4);
	    f5 = f(superior);
	    value = (deltaX/90)*(7*f1 + 32*f2 + 12*f3 + 32*f4 + 7*f5);

	    return value;
	}

	public double integralGrau5(double inferior, double superior)
	{
	    double f1;
	    double f2;
	    double f3;
	    double f4;
	    double f5;
	    double f6;
	    double value;
	    double deltaX;

	    deltaX = (superior - inferior);
	    f1 = f(inferior);
	    f2 = f(inferior + deltaX/5);
	    f3 = f(inferior + deltaX*2/5);
	    f4 = f(inferior + deltaX*3/5);
	    f5 = f(inferior + deltaX*4/5);
	    f6 = f(superior);
	    value = (deltaX/288)*(19*f1 + 75*f2 + 50*f3 + 50*f4 + 75*f5 + 19*f6);

	    return value;
	}
	
	public double integralPartitions(double inferior, double superior, double epsilon, int degree)
	{
	    double current_value = 0.1;
	    double last_value = 1000;
	    double number_of_partitions = 1;
	    double deltaX;
	    double sum;
	    double x;
	    double xi;
	    double xp;

	    while(Math.abs((current_value - last_value)/current_value) > epsilon)
	    {
	        number_of_partitions = 2*number_of_partitions;
	        deltaX = (superior - inferior)/number_of_partitions;
	        sum = 0;

	        for (int i=0; i<number_of_partitions; i++)
	        {
	            x = inferior + i*deltaX;
	            xi = x + deltaX;

				switch (degree) {
				case 0:
					// Caso da integral por retângulo!
					xp = inferior + deltaX * i;
					sum = sum + deltaX * f(xp);

					break;

				case 1:
					// Caso da integral por Trapézio
					xi = inferior + i * deltaX;
					sum = sum + deltaX * (f(xi) + f(xi + deltaX)) / 2;

					break;

				case 2:
					// Newton-Cotes com grau 2
					sum = sum + integralGrau2(x, xi);

					break;

				case 3:
					// Newton-Cotes com grau 3
					sum = sum + integralGrau3(x, xi);
					break;

				case 4:
					// Newton-Cotes com grau 4
					sum = sum + integralGrau4(x, xi);
					break;

				case 5:
					// Newton-Cotes com grau 5
					sum = sum + integralGrau5(x, xi);
					break;
				}

			}

	        last_value = current_value;
	        current_value = sum;
	    }
	    return current_value;
	}

}

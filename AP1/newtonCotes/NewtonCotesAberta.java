//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
public class NewtonCotesAberta {

	public static void main(String[] args) {
		NewtonCotesAberta nca = new NewtonCotesAberta();						//lim inf, lim sup, erro, NCA_grau
		System.out.println("\nO valor da integral é: " + nca.integralPartitions(0, 2, 0.001, 3) + "\n"); 
	}
	
	private double f(double x) {
		// Usando para os testes a função 30x⁴ - 2x³ + 8x - 12.5.
		// A integral dessa função no intervalo de -1 a 1 é -13.
		return (30 * Math.pow(x, 4) - 2 * Math.pow(x, 3) + 8 * x - 12.5);
	}

	public double integralGrau1(double inferior, double superior) {
		double f0;
		double f1;
		double f2;
		double f3;
		double value;
		double deltaX;

		deltaX = (superior - inferior);
		f0 = f(inferior);
		f1 = f(inferior + deltaX / 3);
		f2 = f(inferior + deltaX * 2 / 3);
		f3 = f(superior);
		value = (deltaX / 2) * (0 * f0 + f1 + f2 + 0 * f3);

		return value;
	}

	public double integralGrau2(double inferior, double superior) {
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double value;
		double deltaX;

		deltaX = (superior - inferior);
		f0 = f(inferior);
		f1 = f(inferior + deltaX / 4);
		f2 = f(inferior + deltaX / 2);
		f3 = f(inferior + deltaX * 3 / 4);
		f4 = f(superior);
		value = (deltaX / 3) * (0 * f0 + 2 * f1 - 1 * f2 + 2 * f3 + 0 * f4);

		return value;
	}

	public double integralGrau3(double inferior, double superior) {
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double f5;
		double value;
		double deltaX;

		deltaX = (superior - inferior);
		f0 = f(inferior);
		f1 = f(inferior + deltaX / 5);
		f2 = f(inferior + deltaX * 2 / 5);
		f3 = f(inferior + deltaX * 3 / 5);
		f4 = f(inferior + deltaX * 4 / 5);
		f5 = f(superior);
		value = (deltaX / 24) * (0 * f0 + 11 * f1 + 1 * f2 + 1 * f3 + 11 * f4 + 0 * f5);

		return value;
	}

	public double integralGrau4(double inferior, double superior) {
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double f5;
		double f6;
		double value;
		double deltaX;

		deltaX = (superior - inferior);
		f0 = f(inferior);
		f1 = f(inferior + deltaX / 6);
		f2 = f(inferior + deltaX / 3);
		f3 = f(inferior + deltaX / 2);
		f4 = f(inferior + deltaX * 2 / 3);
		f5 = f(inferior + deltaX * 5 / 6);
		f6 = f(superior);
		value = (deltaX / 20) * (0 * f0 + 11 * f1 - 14 * f2 + 26 * f3 - 14 * f4 + 11 * f5 + 0 * f6);

		return value;
	}

	public double integralGrau5(double inferior, double superior) {
		double f0;
		double f1;
		double f2;
		double f3;
		double f4;
		double f5;
		double f6;
		double f7;
		double value;
		double deltaX;

		deltaX = (superior - inferior);
		f0 = f(inferior);
		f1 = f(inferior + deltaX / 7);
		f2 = f(inferior + deltaX * 2 / 7);
		f3 = f(inferior + deltaX * 3 / 7);
		f4 = f(inferior + deltaX * 4 / 7);
		f5 = f(inferior + deltaX * 5 / 7);
		f6 = f(inferior + deltaX * 6 / 7);
		f7 = f(superior);
		value = (deltaX / 1440) * (0 * f0 + 611 * f1 - 453 * f2 + 562 * f3 + 562 * f4 - 453 * f5 + 611 * f6 + 0 * f7);

		return value;
	}

	public double integralPartitions(double inferior, double superior, double epsilon, int degree) {
		double current_value = 0.1;
		double last_value = 1000;
		double number_of_partitions = 1;
		double deltaX;
		double sum;
		double x;
		double xi;

		while (Math.abs((current_value - last_value) / current_value) > epsilon) {
			number_of_partitions = 2 * number_of_partitions;
			deltaX = (superior - inferior) / number_of_partitions;
			sum = 0;

			for (int i = 0; i < number_of_partitions; i++) {
				x = inferior + i * deltaX;
				xi = x + deltaX;

				switch (degree) {
				case 1:
					// Caso da integral por Trapézio
					sum = sum + integralGrau1(x, xi);

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

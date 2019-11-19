from PIL import Image

#Funcao que decide qual filosofia de derivada utilizar no respectivo pixel.
def derivar(pixel, image_size):		
	
	width, height = image_size
	final_pixel = Image.new('RGB',(width,height), color=0)

	#Considerando a variacaoo de x = 1
	for i in range(width):
		for j in range(height):
			
			if i == 0: #Forward
				temp_pixel = pixel[i+1,j]-pixel[i,j]

			elif i == width-1: #Backward
				temp_pixel = pixel[i,j]-pixel[i-1,j]
				

			else : #Central
				temp_pixel = (pixel[i+1,j]-pixel[i-1,j])/2
	
			final_pixel.putpixel((i,j),int(temp_pixel))
	
	return final_pixel

#Transformando a imagem original em GrayScale.
imagem = Image.open('crash.jpg')
imagem.show()
imagem = imagem.convert('1',None)	
imagem.save('crashGray.jpg')
imagem.show()

#Atribuindo valor aos pixels com load().
pixel = imagem.load()

#Aplicando as derivadas na imagem convertida.
imagem = derivar(pixel,imagem.size)

#Salvando a imagem final e abrindo visualizador de imagens.
imagem.save('crashFinal.jpg')
imagem.show()

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define R 50  //Aqui temos o numero de vezes que faremos GAUSS-SEIDEL

struct CELULA//Nessa linha definimos a celula
{
	int dado;
	int coluna;
	struct CELULA *prox;
};


struct MATRIZ //Nessa linha definimos a matriz
{
	int linha;
	struct CELULA *line;
	struct MATRIZ *prox;
};


int remover(struct MATRIZ *raiz, int lin, int col)//Aqui removemos um valor de poisicao lin x col
{
	struct CELULA *a = NULL;
	struct CELULA *b = NULL;
	struct CELULA *c = NULL;
							//criamos auxiliares para identificar o lugar da matriz 	
	struct MATRIZ *pont = raiz;
	while(pont->linha != lin)
	{
		pont = pont->prox;//procuramos o lugar certo na matriz
	}
	a = pont->line;
	b = pont;	
	if(a == NULL)
	{
		return 0;
	}
	else
	{
		while(a->coluna > col)//temos que se o auxiliar tiver valor maior que o desejado entao ele passou a posicao do valor a ser removido
		{
			if(a->prox != NULL) 
			{
				b = a;
				a = a->prox;
			}
			else 
			{
				return 0;
			}
		}
		if(a->coluna == col) 
		{
			b->prox = a->prox;
			free(a);
			return 1;
		}
		else 
		{
			return 0;
		}
	}
}


void inserir(int dado, int lin, int col, struct MATRIZ *raiz)//inserimos um novo valor a matriz
{
	struct CELULA *a = NULL;
	struct CELULA *b = NULL;
	struct MATRIZ *pont = raiz;
	int i = 0; 
	while(pont->linha != lin) //Aqui achamos o lugar onde deve ser inserido o novo valor
	{
		pont = pont->prox;
	}
	if (pont->line == NULL) 
	{
		a = (struct CELULA*) malloc(sizeof(struct CELULA)); //alocamos memoria para a nova celula
		a->dado = dado;
		a->coluna = col; 
		a->prox = NULL; 
		pont->line = a;
	}
	else 
	{
		a = pont->line;
		if (a->coluna > col) 
		{
			b = (struct CELULA*) malloc(sizeof(struct CELULA));
			b->dado = dado;
			b->coluna = col;
			b->prox = a;
			pont->line = b;
		}
		else
		{
			while(a->prox != NULL && a->coluna < col) 
			{
				a = a->prox;
			}
			b = (struct CELULA*) malloc(sizeof(struct CELULA));
			b->dado = dado;
			b->coluna = col;
			b->prox = a->prox; 
			a->prox = b;
		}
	}
}


int find_position(struct MATRIZ *raiz, int m, int n)
{
	struct MATRIZ *pont = raiz;
	struct CELULA *a = NULL;
	while(pont->linha != m) 
	{
		pont = pont->prox;
	}
	if(pont->line == NULL) 
	{
		return 0; 
	}
	a = pont->line;
	while(a->coluna < n) 
	{
		a = a->prox;
	}
	if(a->coluna == n)
	{
		return(a->dado); 
	}
	else if(a->coluna != n)
	{
		return 0;
	}
}


int soma_linha(struct MATRIZ *raiz, int k) //somamos os valores de uma determinada linha k
{
	struct MATRIZ *pont = raiz;
	struct CELULA *a = NULL;
	int soma = 0;
	while(pont->linha != k) 
	{
		pont = pont->prox;
	}
	if (pont->line == NULL)
	{
		return 0;
	}
	a = pont->line;
	while(a != NULL) 
	{
		soma = soma + a->dado;
		a = a->prox;
	}
	return (soma); 
}


int soma_coluna(struct MATRIZ *raiz, int k)//somamos os valores de uma determinada coluna k
{
	struct MATRIZ *pont = raiz;
	struct CELULA *a = pont->line;
	int soma = 0;
	while (pont != NULL)
	{
		a = pont->line;
		if(a != NULL)
		{ 
			while(a->coluna < k && a != NULL) 
			{
				a = a->prox;
			}
			if(a->coluna == k) 
			{
				soma = soma + a->dado;
			}
		}
		pont = pont->prox;
	}
	return (soma);
}


int diag(struct MATRIZ *raiz, int n)//vamos se a matriz eh diagonal dominante
{
	struct MATRIZ *pont = raiz;
	struct CELULA *a = NULL;
	int c = 0;
	while (pont != NULL)
	{
		a = pont->line;
		while(a != NULL) 
		{
			if(pont->linha == a->coluna)
			{
				c = c - abs(a->dado);
			}
			else if(pont->linha != a->coluna && a->coluna != n-1) 
			{
				c = c + abs(a->dado);
			}
			a = a->prox;
		}
		if(c > 0) 
		{
			printf("\nEssa matriz nao é diagonal dominante, logo não há como utilizar o metodo de convergencia dos metodos iterativos.\n");
			return 0;
		}
		pont = pont->prox;
	}
	if(c <= 0)
	{
		return(1);
	}
}


float aprox(struct MATRIZ *raiz, float *vet, int k, int m, int n)//calculamos os valores aproximados de cada incognita pelo metodo de gauss
{
	float c = find_position(raiz, k, n-1); 
	float b = 0;
	int i;
	struct CELULA *a = raiz->line;
	for(i = 0; i<n -1; i++) 
	{
		if(a->coluna != k && i == a->coluna)
		{
			b = b + (a->dado)*vet[i];
		}
		a = a->prox; 
	}
	b = c - b; 
	b = b/find_position(raiz, k, k);
	return(b);
}


void gauss(struct MATRIZ *raiz, int m, int n)
{
	int i, j;
	float vet[m]; 
	struct MATRIZ *pont = raiz;
	for(i = 0; i < m; i++)
	{
		vet[i] = 0;
	}
	for(j = 0; j <= R; j++) 
	{
		pont = raiz;
		for(i = 0; i < m; i++)
		{
			vet[i] = aprox(pont, vet, i, m, n); 
			pont = pont->prox; 
		}
	}
	printf("\nSolucao:\n");
	for(i = 0; i < m; i++) 
	{
		printf("\nX%d: %.3f", i, vet[i]);
	}
}


int main()//main
{
	int a, b, c, m, n, i, j, opcao, quantidade, count = 0;
	printf("Este programa serve para fazer uma matriz esparsa\n");
	printf("Quantas linhas e quantas colunas terá essa matriz?\n");
	scanf("%d%d", &m, &n);

	struct MATRIZ *raiz = (struct MATRIZ *)malloc(sizeof(struct MATRIZ));
	struct MATRIZ *pont = raiz;
	raiz->prox = NULL;
	raiz->line = NULL;
	raiz->linha = 0;
	for (i = 1; i < m; i++)
	{
		pont->prox = (struct MATRIZ *)malloc(sizeof(struct MATRIZ));
		pont = pont->prox;
		pont->line = NULL;
		pont->prox = NULL;
		pont->linha = i; 
	}
	while (opcao != 7) 
	{
		//Aqui temos o menu
		
		printf("\nMenu: \n");	
		printf("1: Colocar valor na matriz\n2: Resetar Matriz\n3: Valor Posicao\n4: Soma dos valores de cada linha\n5: Soma dos valores de cada coluna\n6:Gauss-Seidel\n7 Sair \nOpcao: ");
		scanf("%d", &opcao);
		if (opcao != 1 && count == 0)
		{
			printf("\nNao existe nenhum dado na sua matriz.\nopcao a opcao 1 e aperte 'Enter'");
		}
		else if (opcao < 1 || opcao > 6)
		{
			printf("Favor entrar com uma das opcoes do menu");
		}
		else
		{
			switch(opcao) 
			{
				case 1:
					printf("\nQuantos dados voce quer inserir na matriz?\n");
					scanf("%d", &quantidade);
					
					count++; 
					printf("\nQual valor? Linha? Coluna?\n");
					
					for (j = 1; j <= quantidade; j++)
					{
						printf("\nValor %d: ", j);
						scanf("%d%d%d", &a, &b, &c);
						inserir(a, b, c, raiz);
					}
					break;
				
				case 2:
					printf("Qual linha? Qual coluna?\n");
					scanf("%d%d", &a, &b);
					c = remover(raiz, a, b);
					break;
				
				case 3:
					printf("Digite os valores da linha e coluna: \n");
					scanf("%d%d", &a, &b);

					printf("\nO dado na posicao %dx%d eh %d.", a, b, find_position(raiz, a, b));
					break;
				
				case 4:
					printf("\nDigite a linha que deseja somar.\n");
					scanf("%d", &a);
					printf("\nA soma dos dados da linha %d eh %d", a, soma_linha(raiz, a));
					break;
				
				case 5:
					printf("\nDigite coluna que deseja somar.\n");
					scanf("%d", &b);
					printf("\nA soma dos dados da coluna %d eh %d", b, soma_coluna(raiz, b));
					break;
					
				case 6:
					if(n != m+1) 
					{
						printf("A matriz nao tem formato n x n+1");
						break;
					}
					if(diag(raiz, n) == 0)
					{
						break;
					}
					gauss(raiz, m, n);
					break;
			}
		}
	}
	printf("\nAu revoir mon ami");//Hasta la vista baby
}

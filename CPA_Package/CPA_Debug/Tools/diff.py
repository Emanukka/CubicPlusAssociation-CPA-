from typing import Callable
import numpy as np


def grad_1(f: Callable, pos:int) -> Callable:

  def diff_f(args:list) ->float:


      d = 1e-5



      # print('args',args)
      x = args[pos]

      # print(x)

      fwd_args = args*1

      rtd_args = args *1


      # print(fwd_args,rtd_args)

      


      # fwd_args_23[pos] = x + (1j**(2/3))*d
      # fwd_args_83[pos] = x + (1j**(8/3))*d
      fwd_args[pos] = x + 1*d

      # f_at_x = f(*args)


      rtd_args[pos] = x - 1*d


      # f_menos = f(*rtd_args)

      # f_mais_23 = f(*fwd_args_23)
      # f_mais_83 = f(*fwd_args_83)

      f_mais = f(*fwd_args)
      f_menos = f(*rtd_args)

      # print(f_mais,'fwd')

      # return np.imag((f_mais_23 - f_mais_83)/(np.sqrt(3)*d))

      return ((f_mais - f_menos)/(2*d))


  return diff_f

def grad_2(f: Callable, pos:int):

    # """
    # f: função a ser diferenciada (ex: np.sin)
    # args: lista de argumentos da função f
    # pos: índice do argumento em relação ao qual estamos diferenciando
    # """

  def diff_f(*args,h=None):

      """
      args: argumentos da função f
      # É passado como key value ( equivale a passar como uma tupla )
      h: passo de discretização
      """
      h = 1e-4
      # Criar cópias dos argumentos para modificar apenas o que precisa
      fwd_args_2 = list(args)  # Para x + 2*h
      fwd_args_1 = list(args)  # Para x + h
      bwd_args_1 = list(args)  # Para x - h
      bwd_args_2 = list(args)  # Para x - 2*h

      # Atualizando o argumento pos com os incrementos e decrementos
      fwd_args_2[pos] = args[pos] + 2 * h
      fwd_args_1[pos] = args[pos] + h
      bwd_args_1[pos] = args[pos] - h
      bwd_args_2[pos] = args[pos] - 2 * h

      # Aproximação da segunda derivada com ordem de precisão 4
      deriv = (-f(*fwd_args_2) + 16*f(*fwd_args_1) - 30*f(*args) +
              16*f(*bwd_args_1) - f(*bwd_args_2)) / (12 * h**2)

      return deriv

  return diff_f


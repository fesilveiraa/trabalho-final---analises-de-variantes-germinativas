# Trabalho-final---analises-de-variantes-germinativas
Trabalho final da pós-graduação da matéria de "Bioinformática aplicada a genômica médica - análises de variantes germinativas e somática"  do Centro de Ensino e Pesquisa Albert Einstein. 

---

Para iniciar o código iremos realizar a integração do Google Drive com o ambiente do Google Colab, permitindo o acesso direto a arquivos armazenados na conta do usuário.

Importando o módulo drive, responsável por gerenciar a conexão entre o Colab e o Google Drive.E montando o Google Drive no diretório /content/drive, tornando seus arquivos acessíveis como se estivessem no sistema de arquivos local do Colab.
O parâmetro force_remount=True força a remontagem do Drive, garantindo que não haja conflitos com conexões anteriores.

```python
from google.colab import drive
import os

# Desmonta o drive se já estiver montado e limpa o diretório
if os.path.isdir('/content/drive'):
    try:
        drive.flush_and_unmount()
    except ValueError:
        pass # Drive not mounted, ignore error
    os.system('rm -rf /content/drive/*') # Remove any residual files

# Cria o diretório de montagem se não existir e o esvazia
if not os.path.exists('/content/drive'):
    os.makedirs('/content/drive')
else:
    os.system('rm -rf /content/drive/*')

drive.mount('/content/drive', force_remount=True)
```
e o devido output deve ser:

```
Mounted at /content/drive
```



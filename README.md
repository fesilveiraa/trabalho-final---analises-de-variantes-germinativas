# Trabalho-final---analises-de-variantes-germinativas
Trabalho final da pós-graduação da matéria de "Bioinformática aplicada a genômica médica - análises de variantes germinativas e somática"  do Centro de Ensino e Pesquisa Albert Einstein. 

---

Inicialmente, é realizada a integração do Google Drive com o ambiente do Google Colab, possibilitando o acesso direto aos arquivos armazenados na conta do usuário. Para isso, utiliza-se o módulo drive, responsável por estabelecer a conexão entre o Colab e o Google Drive. Em seguida, o Drive é montado no diretório /content/drive, permitindo que os arquivos sejam acessados como se estivessem no sistema de arquivos local do Colab. O parâmetro force_remount=True é utilizado para forçar a remontagem do Drive, assegurando que não ocorram conflitos com conexões previamente estabelecidas.

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


---

Este próximo bloco de código é utilizado para atualizar o sistema de pacotes do ambiente Colab e instalar ferramentas essenciais de bioinformática necessárias para a execução do pipeline de análise.

apt-get install -y bwa

Instala o BWA (Burrows-Wheeler Aligner), uma ferramenta amplamente utilizada para o alinhamento de sequências de DNA contra um genoma de referência.

apt-get install -y samtools

Instala o SAMtools, utilizado para manipulação, ordenação, indexação e análise de arquivos SAM/BAM, que armazenam dados de alinhamento.

apt-get install -y bedtools

Instala o BEDTools, um conjunto de utilitários para operações com intervalos genômicos, como interseções, comparações e extração de regiões específicas do genoma.

```bash
%%bash

apt-get update

apt-get install -y bwa
apt-get install -y samtools
apt-get install -y bedtools

echo #Instalando ferramentas complementares..."
apt-get update

apt-get install -y bcftools
apt-get install -y vcftools
```

Baixar a versão 4.1.8.1 do GATK a partir do repositório oficial, extrai seus arquivos para uso no ambiente Colab e remove o arquivo compactado para otimizar o espaço em disco. O GATK é utilizado em análises genômicas e pipelines de chamada de variantes.

```bash
%%bash

wget -q https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip
unzip -q gatk-4.1.8.1.zip
rm gatk-4.1.8.1.zip
```

Executar comandos em ambiente Bash e realiza o download do Picard Tools (versão 2.24.2) diretamente do repositório oficial do Broad Institute. O Picard é utilizado para manipulação e validação de arquivos BAM.

```bash
%%bash

wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar
```

Definir o diretório principal do projeto no Google Drive e cria a estrutura de pastas para armazenar genomas de referência (hg19 e hg38) e arquivos FASTQ. O uso da opção -p garante que as pastas sejam criadas sem erro caso já existam, facilitando a organização e reprodutibilidade do pipeline.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

mkdir -p $MeuDrive/referencia/hg38
mkdir -p $MeuDrive/referencia/hg19
mkdir -p $MeuDrive/dados/fastq
```

Realiza o download do arquivo compactado do cromossomo 10 do genoma humano hg19 a partir do UCSC Genome Browser, descompacta o conteúdo em tempo real e salva o arquivo FASTA no diretório de referência do projeto. Essa sequência é utilizada como genoma de referência nas etapas subsequentes de alinhamento e análise.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"


curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr10.fa.gz" | \
   gunzip -c > "$MeuDrive/referencia/hg19/hg19.fasta"
```
Exibe as 10 primeiras linhas do arquivo FASTA do genoma de referência hg19, permitindo verificar se o download e a descompactação foram realizados corretamente, bem como confirmar o formato do arquivo antes de sua utilização no pipeline.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

head -n 10 $MeuDrive/referencia/hg19/hg19.fasta
```
E o output deve ser similar a esse:

```
>chr10
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```








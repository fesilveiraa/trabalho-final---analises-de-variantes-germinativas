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
Cria os arquivos (5 arquivos auxiliares .amb, .ann, .bwt, .pac, .sa) de índice do BWA a partir do genoma de referência hg19 (cromossomo 10). A opção -a bwtsw define o algoritmo de indexação recomendado para sequências genômicas longas (>2GB), permitindo o alinhamento eficiente das leituras de sequenciamento.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

bwa index \
  -a bwtsw \
  $MeuDrive/referencia/hg19/hg19.fasta
```

Gera o índice FASTA (.fai) do genoma de referência hg19 utilizando o SAMtools. Esse índice permite acesso rápido a regiões específicas do genoma, sendo essencial para etapas posteriores de alinhamento, visualização e análise dos dados de sequenciamento (visualizadores como IVG).

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

samtools faidx $MeuDrive/referencia/hg19/hg19.fasta
```

Remove versões anteriores do arquivo de dicionário, caso existam, e cria o arquivo .dict do genoma de referência hg19 utilizando o Picard. Esse dicionário contém informações estruturais das sequências (nomes, tamanhos e ordem dos cromossomos) e é obrigatório para ferramentas como o GATK.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

rm -f $MeuDrive/referencia/hg19/hg19.dict

java -jar picard.jar CreateSequenceDictionary \
REFERENCE=$MeuDrive/referencia/hg19/hg19.fasta \
OUTPUT=$MeuDrive/referencia/hg19/hg19.dict
```

Este script realiza uma verificação automática da preparação do genoma de referência hg19, conferindo a presença de todos os arquivos essenciais gerados nas etapas anteriores: arquivo FASTA, índice do SAMtools (.fai), dicionário do Picard (.dict) e os arquivos de índice do BWA. Ao final, é exibido um resumo do status da preparação, indicando se o genoma está completo e pronto para ser utilizado nas análises de alinhamento e chamada de variantes.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🔍 Verificação final da preparação do genoma:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

arquivos_essenciais=(
    "hg19.fasta:Genoma FASTA"
    "hg19.fasta.fai:Índice samtools"
    "hg19.dict:Dicionário Picard"
    "hg19.fasta.amb:BWA .amb"
    "hg19.fasta.ann:BWA .ann"
    "hg19.fasta.bwt:BWA .bwt"
    "hg19.fasta.pac:BWA .pac"
    "hg19.fasta.sa:BWA .sa"
)

total=0
presentes=0

for item in "${arquivos_essenciais[@]}"; do
    arquivo=$(echo $item | cut -d: -f1)
    descricao=$(echo $item | cut -d: -f2)

    if [ -f "$MeuDrive/referencia/hg19/$arquivo" ]; then
        tamanho=$(du -h "$MeuDrive/referencia/hg19/$arquivo" | cut -f1)
        echo "✅ $descricao ($tamanho)"
        ((presentes++))
    else
        echo "❌ $descricao - AUSENTE"
    fi
    ((total++))
done

echo ""
echo "📊 RESUMO: $presentes/$total arquivos presentes"

if [ $presentes -eq $total ]; then
    echo "🎉 PREPARAÇÃO COMPLETA! Genoma pronto para uso."
else
    echo "⚠️ Alguns arquivos estão faltando. Revise as etapas anteriores."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

O output ideal deve ser similar a esse: 

```
🔍 Verificação final da preparação do genoma:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✅ Genoma FASTA (132M)
✅ Índice samtools (512)
✅ Dicionário Picard (512)
✅ BWA .amb (512)
✅ BWA .ann (512)
✅ BWA .bwt (130M)
✅ BWA .pac (33M)
✅ BWA .sa (65M)

📊 RESUMO: 8/8 arquivos presentes
🎉 PREPARAÇÃO COMPLETA! Genoma pronto para uso.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

Executa um teste funcional da indexação do genoma hg19, extraindo uma região específica do cromossomo 10 utilizando o samtools faidx. A exibição correta da sequência confirma que o arquivo FASTA e seu índice foram gerados adequadamente e que o genoma está pronto para uso nas etapas de alinhamento e análise do pipeline.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🧪 Teste: Extraindo região chr10:1000-1100"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

samtools faidx "$MeuDrive/referencia/hg19/hg19.fasta" chr10:57227864-57227930

echo ""
echo "✅ Se você viu a sequência acima, a indexação funcionou perfeitamente!"
echo "📏 Essa região tem exatamente 101 bases (1100-1000+1)"
```

Output:
```
🧪 Teste: Extraindo região chr10:1000-1100
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
>chr10:57227864-57227930
GTGTTAGTAAACAAGCAGTTTCTCAAGAGCAGGGGGGAAAAGTTAGTGACAGAAATATGT
TCAAACA

✅ Se você viu a sequência acima, a indexação funcionou perfeitamente!
📏 Essa região tem exatamente 101 bases (1100-1000+1)
```

Este script verifica a existência do diretório de dados FASTQ e lista os arquivos de sequenciamento presentes. Caso nenhum arquivo seja encontrado, mensagens informativas são exibidas para auxiliar na identificação de problemas de organização ou cópia dos dados. Essa etapa garante que os dados de entrada necessários para o alinhamento estejam disponíveis antes da execução do pipeline.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "📁 Verificando arquivos FASTQ..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -d "$MeuDrive/dados/fastq" ]; then
    echo "📄 Arquivos encontrados:"
    ls -lh "$MeuDrive/dados/fastq/"*.fastq 2>/dev/null || {
        echo "❌ Nenhum arquivo FASTQ encontrado!"
        echo "📝 Verifique se os arquivos foram copiados corretamente."
    }
else
    echo "❌ Diretório dados/fastq não encontrado!"
    echo "📝 Crie a estrutura de diretórios e copie os arquivos FASTQ."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Este script localiza um arquivo FASTQ de leituras reversas (R2), exibe exemplos iniciais das sequências e calcula estatísticas básicas, como número total de reads e tamanho do arquivo. Essa etapa permite verificar a integridade e o formato dos dados de sequenciamento, assegurando que os arquivos de entrada estão adequados para as etapas de alinhamento.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "📄 Analisando arquivo FASTQ R2..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

arquivo_r2=$(find "$MeuDrive/dados/fastq/" -name "*R2*.fastq" | head -1)

if [ -f "$arquivo_r2" ]; then
    echo "📁 Arquivo: $(basename "$arquivo_r2")"
    echo ""
    echo "🔍 Primeiras 12 linhas (3 reads completos):"
    head -12 "$arquivo_r2"

    echo ""
    echo "📊 Estatísticas do arquivo:"
    total_linhas=$(wc -l < "$arquivo_r2")
    total_reads=$((total_linhas / 4))
    echo "• Total de linhas: $(printf "%'d" $total_linhas)"
    echo "• Total de reads: $(printf "%'d" $total_reads)"
    echo "• Tamanho do arquivo: $(du -h "$arquivo_r2" | cut -f1)"
else
    echo "❌ Arquivo FASTQ R2 não encontrado!"
    echo "📝 Verifique se os arquivos estão na pasta correta."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Conta o número total de linhas do arquivo FASTQ de leituras diretas (R1). Como cada read FASTQ é composta por quatro linhas, esse valor permite estimar o número total de leituras de sequenciamento, auxiliando na verificação da integridade e do volume dos dados de entrada.

```bash

%%bash
wc -l /content/drive/MyDrive/TRABALHO_FINAL/dados/fastq/cap-ngse-b-2019-chr10_S1_L001_R1_001.fastq
```

Cria o diretório de saída para arquivos BAM e executa o alinhamento das leituras pareadas (R1 e R2) contra o genoma de referência hg19 utilizando o BWA-MEM. São definidos metadados do grupo de leitura (Read Group), como identificador da amostra, biblioteca e plataforma de sequenciamento, garantindo compatibilidade com ferramentas posteriores (ex.: GATK). O resultado do alinhamento é salvo no formato SAM, contendo as leituras alinhadas ao genoma de referência.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

mkdir -p $MeuDrive/dados/bam

SAMPLE="cap-ngse-b-2019"
Biblioteca="Exoma"
Plataforma="Illumina"

arquivo_r1="$MeuDrive/dados/fastq/cap-ngse-b-2019-chr10_S1_L001_R1_001.fastq"
arquivo_r2="$MeuDrive/dados/fastq/cap-ngse-b-2019-chr10_S1_L001_R2_001.fastq"

bwa mem -K 100000000 \
    -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tLB:$Biblioteca\tPL:$Plataforma" \
    "$MeuDrive/referencia/hg19/hg19.fasta" \
    "$arquivo_r1" \
    "$arquivo_r2" > "$MeuDrive/dados/bam/$SAMPLE.sam"
```

Converte o arquivo SAM gerado no alinhamento para o formato BAM, realiza a ordenação das leituras por coordenada genômica e cria o arquivo de índice (.bai) utilizando o SAMtools. Essas etapas são essenciais para otimizar o desempenho e permitir o uso do arquivo BAM em análises subsequentes, como visualização, chamadas de variantes e processamento com o GATK.

``bash 
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

samtools sort -O bam -o "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" "$MeuDrive/dados/bam/$SAMPLE.sam"
samtools index "$MeuDrive/dados/bam/$SAMPLE.sorted.bam"
```

Exibe as cinco primeiras linhas do arquivo SAM gerado no alinhamento, permitindo verificar o cabeçalho e o formato dos registros, incluindo informações de referência e grupos de leitura. Essa inspeção inicial confirma que o alinhamento foi executado corretamente antes das etapas de processamento do BAM.

``` bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

head -5 "$MeuDrive/dados/bam/$SAMPLE.sam"
```

Output: 
```
@SQ	SN:chr10	LN:135534747
@RG	ID:cap-ngse-b-2019	SM:cap-ngse-b-2019	LB:Exoma	PL:Illumina
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -K 100000000 -R @RG\tID:cap-ngse-b-2019\tSM:cap-ngse-b-2019\tLB:Exoma\tPL:Illumina /content/drive/MyDrive/TRABALHO_FINAL/referencia/hg19/hg19.fasta /content/drive/MyDrive/TRABALHO_FINAL/dados/fastq/cap-ngse-b-2019-chr10_S1_L001_R1_001.fastq /content/drive/MyDrive/TRABALHO_FINAL/dados/fastq/cap-ngse-b-2019-chr10_S1_L001_R2_001.fastq
NB551003:113:HG73YBGX7:1:11101:1071:13957	99	chr10	99237520	60	151M	=	99237651	282	GAACTTAGAAGGCAATAATTTTTGCTATTGAATCCCAGTTATGTCAAGGGGTAGAGACAGAGGAGAATACCAACATGACTGTTATCCCTCACCTGTACATGCACTTCCTGGAAGATGGCTTTAAGCACAGAAACAGCCAGCCCTGGGGGCA	AAAAAEAEEEEEEEEE6EEEEEEEEEEEEEEEEEEAEEEEEEEEEAEEEEEEEE/EEEEEEEEEEEEEEEEE/</EEE/EE<EEEEE/EE<A6EEEEAEE/EEAEEEEEEEE/EAEEE<EEEEEEA/EEEE/EAEE/<6AA<<A<A/A<</	NM:i:0	MD:Z:151	MC:Z:151M	AS:i:151	XS:i:0	RG:Z:cap-ngse-b-2019
NB551003:113:HG73YBGX7:1:11101:1071:13957	147	chr10	99237651	60	151M	=	99237520	-282	AACAGCCAGCCCTGGGGGCAGGGCCACACACAGGCTCTGGGGGAGAGGAGAAGGTACGTGAATACCGAAGGAATTGCAGCANNNNCTCCAAGACACAACCCTACGACAGGCCTAATTAGCTACTGTAAGAATCACAGCATCCTGGTTGAGG	EAEEEAEEEEEEEA<A<A<EA6AAAEAAAA<A<<EE/<66/EEAEEEEEEA<EE/EE/EEEEEEEEEEE/EEEE/EEE<EE####EEEEEEEEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEE/EEEEEEEEEAAAAA	NM:i:4	MD:Z:81T0G0G0C66	MC:Z:151M	AS:i:143	XS:i:0	RG:Z:cap-ngse-b-2019
```




























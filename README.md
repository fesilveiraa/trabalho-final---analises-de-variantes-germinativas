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

```bash 
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

Tenta exibir as primeiras linhas do arquivo BAM ordenado, o que resulta em saída ilegível por se tratar de um arquivo binário. Este comando é usado apenas como uma verificação rápida da existência do arquivo. Para inspeção adequada do conteúdo, devem ser utilizados comandos como samtools view ou samtools view -H.

```bash 
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

head -2 "$MeuDrive/dados/bam/$SAMPLE.sorted.bam"
```

O output deve ser algo parecido com isso:

```
�BCe���N�0ǫ�r텗}�݇l��9�C0�"!��d��׳�6>��Hf�����9�������3?�8����ңD��ႦR�Y^0�5;���G�MIL�z�t��n��Y�L�lLS���u�q�;&�p~�rG��U
}�ζ�'�L���<ya�������[:�pI�`��:��!b~�d�H�G�X=R�|���##���2�2���m�ć�9�t��<�5J��%/Ҝ�56��b��) c����bgt�S�Įi�8 �i���xgo���PRn�����owAtҭw�~�J�
```

Exibe o cabeçalho do arquivo BAM ordenado utilizando o samtools view -H. O cabeçalho contém informações essenciais, como genoma de referência, contigs, grupos de leitura (Read Groups) e versões das ferramentas, sendo importante para validar a consistência dos metadados antes das análises subsequentes.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

samtools view -H "$MeuDrive/dados/bam/$SAMPLE.sorted.bam"
```

Este script exibe os primeiros alinhamentos do arquivo BAM ordenado, após converter o conteúdo binário para o formato SAM com samtools view. Além disso, apresenta um resumo das colunas do formato SAM/BAM, facilitando a interpretação dos campos de alinhamento. Essa etapa é utilizada para verificar a qualidade e a coerência dos alinhamentos antes de prosseguir.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "🔍 Conteúdo do arquivo BAM (primeiros 10 alinhamentos):"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "📊 Colunas do formato SAM/BAM:"
echo "1. QNAME: Nome do read"
echo "2. FLAG: Informações binárias (paired, mapped, etc)"
echo "3. RNAME: Cromossomo de referência"
echo "4. POS: Posição no cromossomo (1-based)"
echo "5. MAPQ: Qualidade do mapeamento"
echo "6. CIGAR: Descrição do alinhamento"
echo "7. RNEXT: Cromossomo do par (paired-end)"
echo "8. PNEXT: Posição do par"
echo "9. TLEN: Tamanho do fragmento"
echo "10. SEQ: Sequência do read"
echo "11. QUAL: Qualidades da sequência"
echo ""

echo "📄 Primeiros 10 alinhamentos:"
samtools view "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" | head -10
```

Converte o arquivo BAM ordenado para o formato BED utilizando o BEDTools, representando os alinhamentos como intervalos genômicos. Em seguida, realiza a fusão (merge) de regiões sobrepostas e a ordenação dos intervalos resultantes.

Essa etapa é utilizada para resumir e organizar as regiões genômicas cobertas pelas leituras, facilitando análises baseadas em intervalos, como avaliação de cobertura e interseção com regiões alvo.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

bedtools bamtobed -i "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" > "$MeuDrive/dados/bam/$SAMPLE.bed"
bedtools merge -i "$MeuDrive/dados/bam/$SAMPLE.bed" > "$MeuDrive/dados/bam/$SAMPLE.merged.bed"
bedtools sort -i "$MeuDrive/dados/bam/$SAMPLE.merged.bed" > "$MeuDrive/dados/bam/$SAMPLE.sorted.bed"
```

Exibe as dez primeiras linhas do arquivo BED ordenado, permitindo verificar os intervalos genômicos derivados dos alinhamentos (cromossomo, início e fim). Essa inspeção confirma que a conversão do BAM para BED e a organização das regiões foram realizadas corretamente.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

head -10 "$MeuDrive/dados/bam/$SAMPLE.sorted.bed"
```
Output:
```
chr10	63440	63621
chr10	77630	77788
chr10	80033	80214
chr10	87016	87167
chr10	87514	87664
chr10	87760	87911
chr10	90578	90692
chr10	91421	91597
chr10	92745	95586
chr10	95679	95830

```

Calcula a cobertura média de leitura para cada região genômica definida no arquivo BED, utilizando o arquivo BAM ordenado como referência. A opção -mean retorna a profundidade média de cobertura por intervalo, gerando um arquivo BED com informações quantitativas de cobertura, útil para avaliação da qualidade do sequenciamento e uniformidade da cobertura.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

bedtools coverage \
    -a "$MeuDrive/dados/bam/$SAMPLE.sorted.bed" \
    -b "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" \
    -mean > "$MeuDrive/dados/bam/$SAMPLE.coverage.bed"
```

Este script realiza uma análise descritiva da cobertura de sequenciamento, exibindo exemplos iniciais das regiões analisadas e calculando estatísticas básicas, como número total de regiões, cobertura média, máxima e mínima. Esses indicadores são fundamentais para avaliar a qualidade, profundidade e uniformidade da cobertura, auxiliando na interpretação dos resultados e na validação do pipeline de análise.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📊 Análise detalhada da cobertura:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "🔍 Primeiras 10 regiões com cobertura:"
head -10 "$MeuDrive/dados/bam/$SAMPLE.coverage.bed"

echo ""
echo "📈 Estatísticas gerais de cobertura:"

total_regioes=$(wc -l < "$MeuDrive/dados/bam/$SAMPLE.coverage.bed")
cobertura_media=$(awk '{sum += $4; count++} END {printf "%.2f", sum/count}' "$MeuDrive/dados/bam/$SAMPLE.coverage.bed")
cobertura_maxima=$(awk '{if($4 > max) max = $4} END {printf "%.2f", max}' "$MeuDrive/dados/bam/$SAMPLE.coverage.bed")
cobertura_minima=$(awk 'NR==1{min=$4} {if($4 < min) min = $4} END {printf "%.2f", min}' "$MeuDrive/dados/bam/$SAMPLE.coverage.bed")

echo "• Total de regiões: $(printf "%'d" $total_regioes)"
echo "• Cobertura média: ${cobertura_media}x"
echo "• Cobertura máxima: ${cobertura_maxima}x"
echo "• Cobertura mínima: ${cobertura_minima}x"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

```

Output:
```
📊 Análise detalhada da cobertura:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🔍 Primeiras 10 regiões com cobertura:
chr10	63440	63621	1.6574585
chr10	77630	77788	1.9113925
chr10	80033	80214	1.6685083
chr10	87016	87167	1.0000000
chr10	87514	87664	1.0000000
chr10	87760	87911	1.0000000
chr10	90578	90692	2.0000000
chr10	91421	91597	1.2954545
chr10	92745	95586	55.0165443
chr10	95679	95830	1.0000000

📈 Estatísticas gerais de cobertura:
• Total de regiões: 43,861
• Cobertura média: 6.48x
• Cobertura máxima: 2975.75x
• Cobertura mínima: 1.00x

```

Filtra as regiões genômicas com cobertura média maior ou igual a 20× a partir do arquivo de cobertura, utilizando awk. O resultado é um novo arquivo BED contendo apenas regiões com profundidade considerada adequada para análises confiáveis, como detecção de variantes e avaliação de qualidade do sequenciamento.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "🎯 Filtrando regiões com cobertura ≥ 20x..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

awk -F "\t" '$4 >= 20 {print $0}' "$MeuDrive/dados/bam/$SAMPLE.coverage.bed" > "$MeuDrive/dados/bam/$SAMPLE.coverage.20x.bed"


echo ""
echo "🔍 Primeiras 10 regiões com cobertura ≥ 20x:"
head -10 "$MeuDrive/dados/bam/$SAMPLE.coverage.20x.bed"
```

Output:

```
🎯 Filtrando regiões com cobertura ≥ 20x...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

🔍 Primeiras 10 regiões com cobertura ≥ 20x:
chr10	92745	95586	55.0165443
chr10	130495	131010	22.6174755
chr10	266976	267603	22.5087719
chr10	284964	286270	21.8047466
chr10	287667	288438	29.6900139
chr10	292521	293649	26.4592190
chr10	294053	295084	26.2007751
chr10	297959	298591	21.9367085
chr10	322964	323723	32.0830040
chr10	326870	327593	42.0497932

```

---
Visualização com IGV
🔬 IGV: Integrative Genomics Viewer

IGV é uma ferramenta poderosa para visualização de dados genômicos:

📊 Visualização de alinhamentos (BAM/SAM)
🧬 Análise de variantes (VCF)
📈 Dados de cobertura (BED, bedGraph)
🔍 Navegação interativa no genoma
🎯 Preparação para IGV ✅ Arquivos Necessários:

Genoma de referência: hg19/GRCh37
Arquivo BAM ordenado: cap-ngse-a-2019.sorted.bam
Índice BAM: cap-ngse-a-2019.sorted.bam.bai
Arquivo de cobertura: cap-ngse-a-2019.coverage.20x.bed
🔧 Passos no IGV: 1️⃣ Configurar Genoma de Referência

Genomes → Load Genome from Server
Selecionar: Human hg19 (GRCh37)
2️⃣ Carregar Arquivo BAM

File → Load from File
Selecionar: cap-ngse-a-2019.sorted.bam
⚠️ Importante: Arquivo .bai deve estar na mesma pasta
3️⃣ Navegar para Região com Cobertura

Usar coordenadas das regiões com cobertura ≥ 20x
Formato de busca: chr8:posição_inicial-posição_final

Verificação automática de pré-requisitos do pipeline

Este bloco executa um script para validar a presença de todos os arquivos necessários antes da etapa de chamada de variantes. É definido o diretório de trabalho no Google Drive, garantindo consistência nos caminhos utilizados. Em seguida, é criada uma função para verificar a existência de arquivos críticos, retornando também o tamanho de cada arquivo como forma de validação adicional.

O script checa os principais insumos do pipeline, incluindo o genoma de referência (FASTA, índice .fai e dicionário .dict) e os dados de sequenciamento alinhados (arquivo BAM ordenado e seu índice .bai). Ao final, é apresentado um resumo quantitativo dos arquivos encontrados. Caso todos os pré-requisitos estejam presentes, o pipeline é liberado para execução; caso contrário, o usuário é orientado a executar as etapas anteriores de preparação e alinhamento, garantindo integridade e reprodutibilidade da análise.


```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🔍 Verificação completa de pré-requisitos..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Função para verificar arquivo
verificar_arquivo() {
    local arquivo="$1"
    local descricao="$2"

    if [ -f "$arquivo" ]; then
        local tamanho=$(du -h "$arquivo" | cut -f1)
        echo "✅ $descricao ($tamanho)"
        return 0
    else
        echo "❌ $descricao - AUSENTE"
        return 1
    fi
}

echo "📂 1. Genoma de Referência:"
total=0
presentes=0

if verificar_arquivo "$MeuDrive/referencia/hg19/hg19.fasta" "Genoma FASTA"; then ((presentes++)); fi; ((total++))
if verificar_arquivo "$MeuDrive/referencia/hg19/hg19.fasta.fai" "Índice samtools"; then ((presentes++)); fi; ((total++))
if verificar_arquivo "$MeuDrive/referencia/hg19/hg19.dict" "Dicionário Picard"; then ((presentes++)); fi; ((total++))

echo ""
echo "📊 2. Dados de Sequenciamento:"
if verificar_arquivo "$MeuDrive/dados/bam/cap-ngse-b-2019.sorted.bam" "BAM ordenado"; then ((presentes++)); fi; ((total++))
if verificar_arquivo "$MeuDrive/dados/bam/cap-ngse-b-2019.sorted.bam.bai" "Índice BAM"; then ((presentes++)); fi; ((total++))

echo ""
echo "📋 RESUMO:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Arquivos presentes: $presentes/$total"

if [ $presentes -eq $total ]; then
    echo "🎉 TODOS OS PRÉ-REQUISITOS ATENDIDOS!"
    echo "🚀 Pronto para chamada de variantes."
else
    echo "⚠️ Alguns arquivos estão faltando."
    echo "📝 Execute os notebooks das aulas anteriores primeiro."
    echo ""
    echo "🔗 Ordem recomendada:"
    echo "1. Preparação do Genoma de Referência"
    echo "2. Mapeamento e Alinhamento"
    echo "3. Chamada de Variantes (esta aula)"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Output:
```
🔍 Verificação completa de pré-requisitos...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
📂 1. Genoma de Referência:
✅ Genoma FASTA (132M)
✅ Índice samtools (512)
✅ Dicionário Picard (512)

📊 2. Dados de Sequenciamento:
✅ BAM ordenado (112M)
✅ Índice BAM (185K)

📋 RESUMO:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✅ Arquivos presentes: 5/5
🎉 TODOS OS PRÉ-REQUISITOS ATENDIDOS!
🚀 Pronto para chamada de variantes.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

Análise exploratória e controle de qualidade de arquivos BAM

Este bloco executa uma análise preliminar dos dados de sequenciamento alinhados antes da chamada de variantes. São definidos o diretório de trabalho e o identificador da amostra, garantindo padronização dos nomes de arquivos. O script verifica a existência do arquivo BAM ordenado e, caso esteja presente, exibe seu tamanho em disco como uma checagem inicial de integridade.

Em seguida, o comando samtools flagstat é utilizado para gerar estatísticas básicas do alinhamento, incluindo número total de leituras, leituras mapeadas e taxa de alinhamento. Posteriormente, samtools depth calcula a profundidade de cobertura por posição genômica, permitindo visualizar exemplos iniciais de cobertura e estimar a cobertura média do experimento por meio de um cálculo agregado em awk.

Caso o arquivo BAM não seja encontrado, o script interrompe a análise e orienta o usuário a executar previamente a etapa de mapeamento. Essa verificação assegura a qualidade mínima dos dados e reduz erros nas etapas subsequentes.


```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📊 Análise prévia dos dados BAM..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" ]; then
    echo "📄 Arquivo BAM: $SAMPLE.sorted.bam"
    echo "📏 Tamanho: $(du -h "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" | cut -f1)"

    echo ""
    echo "📈 Estatísticas básicas do BAM:"
    samtools flagstat "$MeuDrive/dados/bam/$SAMPLE.sorted.bam"

    echo ""
    echo "🎯 Região de cobertura (primeiras 5 posições):"
    samtools depth "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" | head -5

    echo ""
    echo "📊 Cobertura média aproximada:"
    samtools depth "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" | \
    awk '{sum+=$3; count++} END {printf "%.1fx (baseado em %d posições)\n", sum/count, count}'

else
    echo "❌ Arquivo BAM não encontrado!"
    echo "📝 Execute o notebook de mapeamento primeiro."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Chamada de variantes germinativas com GATK HaplotypeCaller

Este bloco executa a etapa de chamada de variantes utilizando o GATK HaplotypeCaller, a ferramenta padrão para detecção de SNPs e indels em dados de sequenciamento de nova geração. É criado o diretório de saída para arquivos VCF, garantindo organização e evitando erros por ausência de pasta.

O HaplotypeCaller é executado a partir do genoma de referência no formato FASTA e de um arquivo BAM previamente alinhado e ordenado. O parâmetro --min-base-quality-score 20 filtra bases de baixa qualidade, aumentando a confiabilidade das variantes detectadas. Já o parâmetro --standard-min-confidence-threshold-for-calling 30.0 define um limiar mínimo de confiança (Phred-scaled) para que uma variante seja efetivamente chamada.

Como resultado, é gerado um arquivo VCF contendo as variantes germinativas identificadas na amostra, pronto para etapas posteriores de filtragem, anotação e interpretação.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

mkdir -p $MeuDrive/dados/vcf


./gatk-4.1.8.1/gatk HaplotypeCaller \
    -R "$MeuDrive/referencia/hg19/hg19.fasta" \
    -I "$MeuDrive/dados/bam/$SAMPLE.sorted.bam" \
    -O "$MeuDrive/dados/vcf/$SAMPLE.vcf" \
    --min-base-quality-score 20 \
    --standard-min-confidence-threshold-for-calling 30.0
```

Inspeção e validação inicial do arquivo VCF

Este bloco realiza uma análise exploratória do arquivo VCF gerado na etapa de chamada de variantes, com o objetivo de validar sua integridade antes das etapas de filtragem e anotação. O script verifica a existência do arquivo VCF no diretório de saída e, caso esteja presente, exibe informações básicas como caminho, tamanho em disco e número total de linhas.

Em seguida, é avaliada a estrutura interna do VCF por meio da contagem de linhas de cabeçalho (linhas iniciadas por #) e de registros de variantes propriamente ditos. Essa verificação assegura que o arquivo segue o padrão VCF e contém chamadas de variantes válidas. Por fim, são exibidas as primeiras linhas do cabeçalho, permitindo a conferência manual de metadados críticos, como versão do VCF, parâmetros utilizados e definições de campos INFO e FORMAT.

Caso o arquivo não seja encontrado, o pipeline interrompe a análise e orienta o usuário a executar previamente a etapa de chamada de variantes, garantindo a correta sequência do fluxo de trabalho.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📄 Análise do arquivo VCF gerado..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/vcf/$SAMPLE.vcf" ]; then
    echo "📊 Informações básicas do arquivo:"
    echo "• Localização: $MeuDrive/dados/vcf/$SAMPLE.vcf"
    echo "• Tamanho: $(du -h "$MeuDrive/dados/vcf/$SAMPLE.vcf" | cut -f1)"
    echo "• Total de linhas: $(wc -l < "$MeuDrive/dados/vcf/$SAMPLE.vcf")"

    echo ""
    echo "📋 Estrutura do VCF:"
    linhas_header=$(grep -c '^#' "$MeuDrive/dados/vcf/$SAMPLE.vcf")
    linhas_dados=$(grep -c '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf")
    echo "• Linhas de cabeçalho: $linhas_header"
    echo "• Linhas de dados: $linhas_dados"

    echo ""
    echo "🔍 Cabeçalho do VCF (primeiras 20 linhas):"
    head -20 "$MeuDrive/dados/vcf/$SAMPLE.vcf"

else
    echo "❌ Arquivo VCF não encontrado!"
    echo "📝 Execute a célula de chamada de variantes primeiro."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```














































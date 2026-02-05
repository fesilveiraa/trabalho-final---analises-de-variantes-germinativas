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


Estatísticas descritivas e caracterização das variantes chamadas

Este bloco executa uma análise quantitativa e qualitativa do arquivo VCF gerado, com o objetivo de caracterizar o conjunto de variantes antes das etapas de filtragem avançada e interpretação. O script verifica a existência do VCF e contabiliza o número total de variantes chamadas (linhas não comentadas).

Em seguida, as variantes são classificadas de forma simplificada em SNPs e INDELs, com base no comprimento dos alelos de referência e alternativo. Também é avaliada a distribuição das variantes em diferentes limiares de qualidade (campo QUAL), permitindo uma estimativa rápida da confiabilidade das chamadas.

Por fim, são exibidos exemplos das primeiras variantes identificadas, incluindo os principais campos do VCF (cromossomo, posição, alelos, qualidade e informações), facilitando a inspeção manual e a validação do formato. Caso nenhuma variante seja detectada, o script fornece possíveis explicações técnicas, auxiliando no diagnóstico de problemas experimentais ou de parametrização.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📊 Estatísticas detalhadas das variantes..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/vcf/$SAMPLE.vcf" ]; then

    # Contagem geral
    echo "🔢 Contagens gerais:"
    variantes=$(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | wc -l)
    echo "• Variantes chamadas: $(printf "%'d" $variantes)"

    if [ $variantes -gt 0 ]; then
        echo ""
        echo "🧬 Análise dos tipos de variantes:"

        # Identificar SNPs e INDELs
        snps=$(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | \
               awk 'length($4)==1 && length($5)==1' | wc -l)
        indels=$(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | \
                awk 'length($4)!=length($5)' | wc -l)

        echo "• SNPs (Single Nucleotide Polymorphisms): $snps"
        echo "• INDELs (Inserções/Deleções): $indels"

        # Distribuição por qualidade
        echo ""
        echo "📈 Distribuição por qualidade (QUAL):"
        echo "• QUAL ≥ 30: $(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | awk '$6 >= 30' | wc -l)"
        echo "• QUAL ≥ 50: $(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | awk '$6 >= 50' | wc -l)"
        echo "• QUAL ≥ 100: $(grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | awk '$6 >= 100' | wc -l)"

        echo ""
        echo "🎯 Primeiras 5 variantes identificadas:"
        grep '^[^#]' "$MeuDrive/dados/vcf/$SAMPLE.vcf" | head -5 | \
        cut -f1-8 | column -t

    else
        echo "⚠️ Nenhuma variante identificada."
        echo "💡 Isso pode indicar:"
        echo "   • Baixa cobertura na região"
        echo "   • Parâmetros muito restritivos"
        echo "   • Região conservada no cromossomo 8"
    fi

else
    echo "❌ Arquivo VCF não encontrado!"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Output:
```
📊 Estatísticas detalhadas das variantes...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🔢 Contagens gerais:
• Variantes chamadas: 6,965

🧬 Análise dos tipos de variantes:
• SNPs (Single Nucleotide Polymorphisms): 6412
• INDELs (Inserções/Deleções): 553

📈 Distribuição por qualidade (QUAL):
• QUAL ≥ 30: 6965
• QUAL ≥ 50: 5731
• QUAL ≥ 100: 4655

🎯 Primeiras 5 variantes identificadas:
chr10  80119  .  C  G  78.32    .  AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=25.36;SOR=0.693
chr10  80124  .  A  G  78.32    .  AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=28.73;SOR=0.693
chr10  93581  .  G  T  64.64    .  AC=1;AF=0.500;AN=2;BaseQRankSum=3.246;DP=85;ExcessHet=3.0103;FS=13.366;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=0.76;ReadPosRankSum=2.051;SOR=3.549
chr10  93603  .  C  T  1169.64  .  AC=1;AF=0.500;AN=2;BaseQRankSum=-1.611;DP=104;ExcessHet=3.0103;FS=6.880;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=12.71;ReadPosRankSum=-2.915;SOR=0.290
chr10  93616  .  C  T  245.64   .  AC=1;AF=0.500;AN=2;BaseQRankSum=3.009;DP=108;ExcessHet=3.0103;FS=13.311;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=2.27;ReadPosRankSum=-1.581;SOR=2.527
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

```

Filtragem de variantes com base em critérios de qualidade (QUAL)

Este bloco aplica uma filtragem inicial ao arquivo VCF com o objetivo de selecionar apenas variantes de alta confiabilidade. Após verificar a existência do VCF de entrada, é utilizado o bcftools filter para reter exclusivamente variantes com valor de qualidade (QUAL) igual ou superior a 100, um critério conservador comumente empregado para reduzir falsos positivos.

Em seguida, o script compara quantitativamente o número de variantes antes e após a filtragem, permitindo avaliar o impacto do filtro aplicado e o percentual de variantes mantidas. Essa etapa é essencial para verificar se os critérios de qualidade estão excessivamente restritivos ou adequados ao conjunto de dados analisado.

Por fim, caso variantes sejam retidas, são exibidos exemplos das primeiras variantes filtradas, facilitando a inspeção manual. Se nenhuma variante atender aos critérios, o pipeline fornece sugestões técnicas para ajuste dos parâmetros, como redução do limiar de qualidade ou verificação da cobertura de sequenciamento. Essa abordagem garante transparência, rastreabilidade e controle de qualidade na seleção de variantes para análises subsequentes.

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "🔍 Filtragem de variantes de alta qualidade..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/vcf/$SAMPLE.vcf" ]; then

    echo "📋 Aplicando filtros de qualidade:"
    echo "• QUAL ≥ 100 (qualidade da chamada)"
    echo ""

    # Aplicar filtros de qualidade
    bcftools filter -i 'QUAL>=100' "$MeuDrive/dados/vcf/$SAMPLE.vcf" > "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf"

    echo "✅ Filtragem concluída!"
    echo ""

    # Estatísticas antes e depois da filtragem
    echo "📊 Comparação antes/depois da filtragem:"

    variantes_total=$(bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.vcf" | wc -l)
    variantes_filtradas=$(bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | wc -l)

    echo "• Variantes antes da filtragem: $variantes_total"
    echo "• Variantes após filtragem: $variantes_filtradas"

    if [ $variantes_total -gt 0 ]; then
        percentual=$(awk "BEGIN {printf \"%.1f\", ($variantes_filtradas/$variantes_total)*100}")
        echo "• Percentual mantido: $percentual%"
    fi

    # Mostrar variantes filtradas se existirem
    if [ $variantes_filtradas -gt 0 ]; then
        echo ""
        echo "🎯 Variantes de alta qualidade identificadas:"
        bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | head -10 | \
        cut -f1-8 | column -t
    else
        echo ""
        echo "⚠️ Nenhuma variante passou pelos filtros de qualidade."
        echo "💡 Sugestões:"
        echo "   • Reduzir threshold de qualidade (QUAL < 100)"
        echo "   • Verificar cobertura da região"
    fi

else
    echo "❌ Arquivo VCF não encontrado!"
    echo "📝 Execute a chamada de variantes primeiro."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Interpretação do formato VCF e análise detalhada de variantes

Este bloco tem como objetivo explicar a estrutura do formato VCF (Variant Call Format) e exemplificar a interpretação de uma variante filtrada de alta qualidade. O script valida a existência do arquivo VCF filtrado e apresenta a descrição funcional de cada coluna padrão do VCF, facilitando a compreensão dos campos utilizados ao longo do pipeline.

Em seguida, é exibido o cabeçalho das colunas (#CHROM), garantindo a conferência do layout do arquivo e da ordem dos campos. Caso o VCF contenha variantes, o script seleciona a primeira entrada e realiza uma decomposição dos principais atributos, incluindo localização genômica, alelos de referência e alternativo, qualidade da chamada e status do filtro.

A variante é então classificada automaticamente como SNP, inserção ou deleção, com base no comprimento relativo dos alelos REF e ALT. Por fim, a linha completa da variante é exibida, permitindo a correlação direta entre a interpretação didática e o registro bruto do VCF. Essa etapa é fundamental para consolidar o entendimento do formato VCF e preparar o usuário para análises de anotação funcional e interpretação clínica.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📋 Interpretação detalhada do formato VCF..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" ]; then

    echo "📖 Estrutura das colunas VCF:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "1. CHROM  : Cromossomo"
    echo "2. POS    : Posição (1-based)"
    echo "3. ID     : Identificador (rs number)"
    echo "4. REF    : Alelo de referência"
    echo "5. ALT    : Alelo alternativo"
    echo "6. QUAL   : Qualidade da chamada (Phred score)"
    echo "7. FILTER : Status do filtro (PASS/FAIL)"
    echo "8. INFO   : Informações adicionais"
    echo "9. FORMAT : Formato dos dados da amostra"
    echo "10. SAMPLE: Dados específicos da amostra"
    echo ""

    # Mostrar linha de header das colunas
    echo "📄 Cabeçalho das colunas:"
    grep '^#CHROM' "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf"
    echo ""

    # Análise de uma variante específica (se existir)
    if [ $(bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | wc -l) -gt 0 ]; then
        echo "🔍 Análise detalhada da primeira variante:"
        variante=$(bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | head -1)

        # Separar campos
        chrom=$(echo "$variante" | cut -f1)
        pos=$(echo "$variante" | cut -f2)
        ref=$(echo "$variante" | cut -f4)
        alt=$(echo "$variante" | cut -f5)
        qual=$(echo "$variante" | cut -f6)
        filter=$(echo "$variante" | cut -f7)

        echo "• Localização: $chrom:$pos"
        echo "• Mudança: $ref → $alt"
        echo "• Qualidade: $qual"
        echo "• Status: $filter"

        # Determinar tipo de variante
        if [ ${#ref} -eq 1 ] && [ ${#alt} -eq 1 ]; then
            echo "• Tipo: SNP (Single Nucleotide Polymorphism)"
        elif [ ${#ref} -gt ${#alt} ]; then
            echo "• Tipo: Deleção ($(( ${#ref} - ${#alt} )) base(s))"
        elif [ ${#ref} -lt ${#alt} ]; then
            echo "• Tipo: Inserção ($(( ${#alt} - ${#ref} )) base(s))"
        fi

        echo ""
        echo "📊 Linha completa da variante:"
        echo "$variante"

    else
        echo "ℹ️ Nenhuma variante disponível para análise detalhada."
    fi

else
    echo "❌ Arquivo VCF filtrado não encontrado!"
    echo "📝 Execute a filtragem primeiro."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Análise estatística de variantes com VCFtools (razão Ts/Tv)

Este bloco executa uma análise estatística do arquivo VCF filtrado utilizando o VCFtools, com foco no cálculo da razão entre transições (Ts) e transversões (Tv). O script verifica a existência do VCF de entrada para garantir a correta execução da etapa.

O comando --TsTv-summary gera um resumo das substituições nucleotídicas observadas, contabilizando eventos de transição (A↔G, C↔T) e transversão (purina↔pirimidina). A razão Ts/Tv é amplamente utilizada como métrica de controle de qualidade em análises genômicas, pois valores esperados (≈2.0–3.0 para dados humanos germinativos) indicam boa qualidade das chamadas de variantes.

Os resultados são salvos em um arquivo de saída e exibidos no terminal para inspeção imediata. Essa etapa fornece uma validação adicional da confiabilidade do conjunto de variantes filtradas antes das fases de anotação funcional e interpretação biológica ou clínica.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "📈 Análise estatística detalhada com vcftools..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" ]; then

    echo "🔍 Executando análises estatísticas..."

    # Estatísticas gerais
    vcftools --vcf "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" --TsTv-summary --out "$MeuDrive/dados/vcf/$SAMPLE"
    cat  "$MeuDrive/dados/vcf/$SAMPLE.TsTv.summary"

else
    echo "❌ Arquivo VCF filtrado não encontrado!"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

 8. Visualização e Validação 🔬 Preparação para IGV

Para validar as variantes chamadas, podemos visualizá-las no IGV junto com os dados de alinhamento. 📋 Arquivos Necessários para IGV:

Genoma de referência: hg38
Arquivo BAM: cap-ngse-a-2019.sorted.bam + .bai
Arquivo VCF: cap-ngse-a-2019.filtered.vcf
🎯 Instruções para IGV:

Carregar genoma: Human hg38
Carregar BAM: Para ver alinhamentos
Carregar VCF: Para ver variantes
Navegar para variantes: Usar coordenadas das variantes

Preparação e validação de arquivos para visualização no IGV

Este bloco realiza a verificação final dos arquivos necessários para inspeção visual das variantes no IGV, uma etapa fundamental de validação manual. O script checa a presença do arquivo BAM alinhado e ordenado, seu respectivo índice (.bai) e o arquivo VCF filtrado, garantindo compatibilidade com o IGV.

Quando todos os arquivos estão disponíveis, o script identifica automaticamente regiões genômicas contendo variantes de alta qualidade e sugere janelas expandidas ao redor dessas posições, facilitando a navegação e o zoom no IGV para avaliação do suporte por leituras individuais. Essa abordagem reduz o tempo de inspeção manual e direciona o usuário para regiões informativas.

Por fim, são fornecidas instruções práticas para carregamento do genoma de referência, arquivos BAM e VCF no IGV, orientando o fluxo de validação visual. Caso algum arquivo esteja ausente, o pipeline interrompe a etapa e orienta a execução prévia das fases necessárias, assegurando a integridade e a reprodutibilidade da análise.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"
SAMPLE="cap-ngse-b-2019"

echo "👁️ Preparação para visualização no IGV..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Verificar arquivos necessários para IGV
echo "✅ Checklist de arquivos para IGV:"

arquivos_igv=(
    "$MeuDrive/dados/bam/$SAMPLE.sorted.bam:Arquivo BAM"
    "$MeuDrive/dados/bam/$SAMPLE.sorted.bam.bai:Índice BAM"
    "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf:Arquivo VCF"
)

todos_presentes=true
for item in "${arquivos_igv[@]}"; do
    arquivo=$(echo $item | cut -d: -f1)
    descricao=$(echo $item | cut -d: -f2)

    if [ -f "$arquivo" ]; then
        echo "✅ $descricao"
    else
        echo "❌ $descricao - AUSENTE"
        todos_presentes=false
    fi
done

echo ""
if [ "$todos_presentes" = true ]; then
    echo "🎉 Todos os arquivos estão disponíveis para IGV!"

    # Sugerir regiões para visualização
    if [ -f "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" ]; then
        variantes_count=$(bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | wc -l)

        if [ $variantes_count -gt 0 ]; then
            echo ""
            echo "📍 Regiões recomendadas para visualização:"

            bcftools view -H "$MeuDrive/dados/vcf/$SAMPLE.filtered.vcf" | head -5 | \
            while read linha; do
                chrom=$(echo "$linha" | cut -f1)
                pos=$(echo "$linha" | cut -f2)
                ref=$(echo "$linha" | cut -f4)
                alt=$(echo "$linha" | cut -f5)
                qual=$(echo "$linha" | cut -f6)

                # Criar região expandida para visualização
                start=$((pos - 50))
                end=$((pos + 50))

                echo "• $chrom:$start-$end ($ref→$alt, QUAL=$qual)"
            done
        fi
    fi

    echo ""
    echo "🔧 Passos no IGV:"
    echo "1. Genomes → Load Genome from Server → Human hg19"
    echo "2. File → Load from File → Selecionar BAM"
    echo "3. File → Load from File → Selecionar VCF"
    echo "4. Navegar para uma das regiões sugeridas acima"
    echo "5. Zoom in para ver reads individuais"

else
    echo "⚠️ Alguns arquivos estão faltando para visualização no IGV."
    echo "📝 Execute as etapas anteriores do pipeline primeiro."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Output:
```
👁️ Preparação para visualização no IGV...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✅ Checklist de arquivos para IGV:
✅ Arquivo BAM
✅ Índice BAM
✅ Arquivo VCF

🎉 Todos os arquivos estão disponíveis para IGV!

📍 Regiões recomendadas para visualização:
• chr10:93553-93653 (C→T, QUAL=1169.64)
• chr10:93566-93666 (C→T, QUAL=245.64)
• chr10:93644-93744 (T→C, QUAL=367.64)
• chr10:93682-93782 (A→G, QUAL=206.64)
• chr10:93895-93995 (G→A, QUAL=768.64)

🔧 Passos no IGV:
1. Genomes → Load Genome from Server → Human hg19
2. File → Load from File → Selecionar BAM
3. File → Load from File → Selecionar VCF
4. Navegar para uma das regiões sugeridas acima
5. Zoom in para ver reads individuais
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

```

Verificação de pré-requisitos para anotação de variantes

Este bloco realiza uma checagem automatizada dos arquivos VCF necessários para a etapa de anotação de variantes, garantindo que o pipeline seja executado apenas quando os insumos mínimos estiverem disponíveis. É definido o diretório de trabalho e criada uma função reutilizável para verificar a existência dos arquivos, exibindo também o tamanho em disco como validação adicional de integridade.

O script verifica a presença tanto do VCF bruto gerado pela chamada de variantes quanto do VCF previamente filtrado, contabilizando quantos arquivos estão disponíveis. Caso ao menos um VCF esteja presente, o pipeline é liberado para a etapa de anotação, priorizando o uso do arquivo filtrado quando disponível.

Por fim, são exibidas estatísticas básicas do arquivo selecionado, incluindo o número total de variantes a serem anotadas. Caso nenhum VCF seja encontrado, o script orienta o usuário a executar as etapas anteriores do fluxo de trabalho, assegurando ordem lógica, reprodutibilidade e confiabilidade na anotação funcional e clínica das variantes.

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🔍 Verificação de pré-requisitos para anotação..."
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

echo "📊 Arquivos VCF necessários:"
total=0
presentes=0

# Verificar arquivos VCF das aulas anteriores
arquivos_vcf=(
    "$MeuDrive/dados/vcf/cap-ngse-b-2019.filtered.vcf:VCF filtrado"
    "$MeuDrive/dados/vcf/cap-ngse-b-2019.vcf:VCF de variantes"
)

for item in "${arquivos_vcf[@]}"; do
    arquivo=$(echo $item | cut -d: -f1)
    descricao=$(echo $item | cut -d: -f2)

    if verificar_arquivo "$arquivo" "$descricao"; then
        ((presentes++))
    fi
    ((total++))
done

echo ""
echo "📋 RESUMO:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Arquivos presentes: $presentes/$total"

if [ $presentes -gt 0 ]; then
    echo "🎉 Pelo menos um arquivo VCF disponível!"
    echo "🚀 Pronto para anotação de variantes."

    # Mostrar estatísticas básicas do melhor arquivo disponível
    if [ -f "$MeuDrive/dados/vcf/cap-ngse-b-2019.filtered.vcf" ]; then
        arquivo_trabalho="$MeuDrive/dados/vcf/cap-ngse-b-2019.filtered.vcf"
        echo "📄 Arquivo principal: cap-ngse-b-2019.filtered.vcf"
    elif [ -f "$MeuDrive/dados/vcf/cap-ngse-b-2019.vcf" ]; then
        arquivo_trabalho="$MeuDrive/dados/vcf/cap-ngse-b-2019.variants.vcf"
        echo "📄 Arquivo principal: cap-ngse-b-2019.vcf"
    fi

    if [ -n "$arquivo_trabalho" ]; then
        variantes=$(grep -c '^[^#]' "$arquivo_trabalho")
        echo "🧬 Variantes a anotar: $variantes"
    fi

else
    echo "⚠️ Nenhum arquivo VCF encontrado!"
    echo "📝 Execute os notebooks das aulas anteriores primeiro:"
    echo "   1. Preparação do Genoma de Referência"
    echo "   2. Mapeamento e Alinhamento"
    echo "   3. Chamada de Variantes"
    echo "   4. Anotação (esta aula)"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Output:
```
🔍 Verificação de pré-requisitos para anotação...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
📊 Arquivos VCF necessários:
✅ VCF filtrado (877K)
✅ VCF de variantes (1.3M)

📋 RESUMO:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✅ Arquivos presentes: 2/2
🎉 Pelo menos um arquivo VCF disponível!
🚀 Pronto para anotação de variantes.
📄 Arquivo principal: cap-ngse-b-2019.filtered.vcf
🧬 Variantes a anotar: 4655
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

Download, instalação e validação do ANNOVAR

Este bloco automatiza a instalação do ANNOVAR, uma das ferramentas mais utilizadas para anotação funcional e genética de variantes. O script verifica se o diretório annovar já existe no ambiente, evitando reinstalações desnecessárias e garantindo idempotência do pipeline.

O script realiza o download da versão pública mais recente diretamente do site oficial, extrai o conteúdo do arquivo compactado (.tar.gz) e remove o arquivo temporário para otimizar o uso de espaço em disco.

Após a instalação, é feita uma verificação de integridade, confirmando a presença do script principal annotate_variation.pl e listando os scripts Perl disponíveis no diretório. Essa checagem garante que o ambiente esteja corretamente preparado para as próximas etapas de anotação de variantes, reduzindo falhas de execução e facilitando o diagnóstico de problemas.

```bash

%%bash

echo "📥 Baixando e instalando ANNOVAR..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Verificar se ANNOVAR já está instalado
if [ -d "annovar" ]; then
    echo "✅ ANNOVAR já está instalado!"
    echo "📁 Localização: $(pwd)/annovar"
else
    echo "📥 Baixando ANNOVAR..."

    # Download do ANNOVAR (versão pública)
    wget -q http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

    echo "📦 Extraindo ANNOVAR..."
    tar -xzf annovar.latest.tar.gz

    echo "🧹 Limpando arquivo temporário..."
    rm annovar.latest.tar.gz

    echo "✅ ANNOVAR instalado com sucesso!"
fi

# Verificar instalação
if [ -f "annovar/annotate_variation.pl" ]; then
    echo ""
    echo "🔍 Verificando instalação:"
    echo "• Script principal: $(ls -la annovar/annotate_variation.pl | awk '{print $1, $5, $9}')"
    echo "• Scripts disponíveis: $(ls annovar/*.pl | wc -l) arquivos"

    echo ""
    echo "📋 Scripts principais do ANNOVAR:"
    ls annovar/*.pl | while read script; do
        nome=$(basename "$script")
        echo "  • $nome"
    done

else
    echo "❌ Erro na instalação do ANNOVAR!"
    echo "📝 Verifique a conectividade e tente novamente."
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
```

Download do banco de dados refGene para anotação com ANNOVAR Este bloco realiza o download do banco de dados refGene, utilizado pelo ANNOVAR para anotação baseada em genes, a partir dos servidores oficiais do ANNOVAR. A variável db_name é definida para facilitar a reutilização do código com outros bancos de dados, promovendo modularidade do pipeline.

O script annotate_variation.pl é executado com a opção -downdb, indicando que a ação é o download do banco de dados, e -buildver hg38, especificando a versão do genoma de referência humano utilizada na análise. O parâmetro -webfrom annovar define a origem dos dados, enquanto o diretório annovar/humandb é utilizado para armazenar localmente os arquivos baixados.

Essa etapa prepara o ambiente para a anotação funcional das variantes, permitindo a associação das posições genômicas às estruturas gênicas conhecidas (genes, éxons, íntrons e regiões regulatórias) de acordo com a versão do genoma selecionada.

```bash

%%bash

db_name="refGene"

perl annovar/annotate_variation.pl -buildver hg19 -downdb \
    -webfrom annovar "$db_name" annovar/humandb
```

```bash

%%bash

MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🗄️ Baixando bancos de dados essenciais..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

db_name="gnomad_exome"

if perl annovar/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar "$db_name" annovar/humandb 2>/dev/null; then
    echo "✅ $description baixado com sucesso"
else
    echo "⚠️ $description - erro no download"
fi

echo ""
echo "📊 Resumo do banco baixado:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


if [ -f "annovar/humandb/hg19_$db_name.txt" ]; then
    tamanho=$(du -h "annovar/humandb/hg19_$db_name.txt" | cut -f1)
    echo "✅ $db_name ($tamanho)"
    echo "🎉 Banco $db_name pronto para anotação!"
else
    echo "❌ $db_name - não disponível"
    echo "⚠️ Problemas no download"
```

```bash

%%bash

MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🗄️ Baixando bancos de dados essenciais..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

db_name="revel"

if perl annovar/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar "$db_name" annovar/humandb/ 2>/dev/null; then
    echo "✅ $description baixado com sucesso"
else
    echo "⚠️ $description - erro no download"
fi

echo ""
echo "📊 Resumo do banco baixado:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


if [ -f "annovar/humandb/hg19_$db_name.txt" ]; then
    tamanho=$(du -h "annovar/humandb/hg19_$db_name.txt" | cut -f1)
    echo "✅ $db_name ($tamanho)"
    echo "🎉 Banco $db_name pronto para anotação!"
else
    echo "❌ $db_name - não disponível"
    echo "⚠️ Problemas no download"
fi
```

```bash

%%bash

MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🗄️ Baixando bancos de dados essenciais..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

db_name="clinvar_20200316"

if perl annovar/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar "$db_name" annovar/humandb/ 2>/dev/null; then
    echo "✅ $description baixado com sucesso"
else
    echo "⚠️ $description - erro no download"
fi

echo ""
echo "📊 Resumo do banco baixado:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


if [ -f "annovar/humandb/hg19_$db_name.txt" ]; then
    tamanho=$(du -h "annovar/humandb/hg19_$db_name.txt" | cut -f1)
    echo "✅ $db_name ($tamanho)"
    echo "🎉 Banco $db_name pronto para anotação!"
else
    echo "❌ $db_name - não disponível"
    echo "⚠️ Problemas no download"
fi
```

```bash
%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "📄 Preparando arquivo VCF para anotação..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Verificar conteúdo do VCF
echo "🧬 Variantes no arquivo: $MeuDrive/dados/vcf/cap-ngse-b-2019.filtered.vcf"

# Criar diretório de anotação
mkdir -p "$MeuDrive/dados/annotation"

perl annovar/convert2annovar.pl -format vcf4 "$MeuDrive/dados/vcf/cap-ngse-b-2019.filtered.vcf" \
    > "$MeuDrive/dados/annotation/variantes.avinput"

# Verificar conversão
if [ -f "$MeuDrive/dados/annotation/variantes.avinput" ]; then
    linhas_convertidas=$(wc -l < "$MeuDrive/dados/annotation/variantes.avinput")
    echo "✅ Conversão concluída!"
    echo "📄 Arquivo gerado: variantes.avinput"
    echo "📊 Linhas convertidas: $linhas_convertidas"
else
    echo "❌ Erro na conversão para formato ANNOVAR!"
fi
```

```bash

%%bash
MeuDrive="/content/drive/MyDrive/TRABALHO_FINAL"

echo "🧬 Executando anotação básica (gene-based)..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -f "$MeuDrive/dados/annotation/variantes.avinput" ]; then
    echo "❌ Arquivo ANNOVAR input não encontrado!"
    echo "📝 Execute a célula de preparação anterior."
    exit 1
fi

echo "🔄 Executando anotação genômica..."

# Anotação básica com RefSeq
perl annovar/annotate_variation.pl -geneanno -buildver hg19 \
    "$MeuDrive/dados/annotation/variantes.avinput" annovar/humandb/
```






























































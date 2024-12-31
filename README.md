# Análise Geoespacial de Solos Usando R

Este projeto realiza uma análise geoespacial de solos, com foco em identificar padrões espaciais e variabilidade dos atributos analisados. Ele foi desenvolvido utilizando ferramentas estatísticas e espaciais em R.

## Relatório

O relatório detalhado sobre este projeto pode ser acessado [aqui](Relatório%20de%20Análise%20Geoespacial%20de%20Solos%20Usando%20R.pdf).

## Estrutura do Projeto

- **`imovel/`**: Contém os shapefiles do imóvel analisado.
- **`localizacao/`**: Arquivos relacionados ao mapa de localização e pontos de amostragem.
- **`saida/`**: Produtos gerados, como mapas e tabelas.
- **`analise_estatistica_argila_organizado.R`**: Script R para a análise estatística e geoespacial.

## Metodologia

1. **Preparação dos Dados**: Dados físico-químicos de solos foram importados e analisados.
2. **Análise Variográfica**: Avaliação da dependência espacial através de semivariogramas.
3. **Interpolação por Krigagem**: Geração de mapas contínuos de atributos do solo.
4. **Mapa de Variabilidade Espacial**: Representação das diferenças relativas nos dados.

## Ferramentas Utilizadas

- [R](https://www.r-project.org/)
  - Pacotes: `sf`, `raster`, `gstat`, `ggplot2`

## Resultados

Os resultados incluem mapas interpolados, análise estatística e visualizações geoespaciais.

## Licença

Este projeto está licenciado sob a [MIT License](LICENSE).

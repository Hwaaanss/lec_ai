# GSE25066 분석 파이프라인

다음 순서로 진행되는 분석 파이프라인입니다.

- 전처리 및 코호트 정의
- EDA 및 품질 점검
- 생물학적 해석이 가능한 feature engineering
- ML 모델링 및 검증
- 요약 리포트 생성

이 파이프라인은 다음과 같이 구성됩니다.

- `R`: 데이터 파싱, annotation, EDA, DE, enrichment, pathway score 계산
- `Python` (`conda activate bc` 환경): ML 모델링 및 성능 평가

## 구성

- `config/analysis_config.yaml`: 파이프라인 설정 파일
- `config/bc_environment.yml`: Python 환경 정의 파일
- `config/immune_signatures.tsv`: 사용자 정의 immune activation/suppression signature
- `scripts/R/install_r_packages.R`: 필요한 R 패키지를 글로벌에 설치
- `scripts/R/01_prepare_data.R`: GEO matrix 파싱, probe annotation, 코호트 생성
- `scripts/R/02_eda.R`: 코호트 요약, QC plot, PCA, 임상 변수 연관성 검정
- `scripts/R/03_bi_analysis.R`: ssGSEA/GSVA, DE, GO BP/KEGG enrichment, 모델 입력 생성
- `scripts/python/train_models.py`: nested CV, source-transfer 평가, 해석 결과 생성
- `scripts/R/04_compile_summary.R`: 결과를 markdown 요약 리포트로 통합
- `main.py`: 순차 실행용 파이프라인 엔트리포인트

## 환경 준비

### 1. Python 환경

필요한 conda 환경을 생성합니다.

```bash
conda env create -n bc -f config/bc_environment.yml
```

이미 환경이 있으면 업데이트합니다.

```bash
conda env update -n bc -f config/bc_environment.yml --prune
```

### 2. R 패키지

부족한 R 패키지를 글로벌에 설치합니다.

```bash
Rscript scripts/R/install_r_packages.R
```

## 실행

```bash
python main.py
```

필요하면 Python 모델링 단계만 건너뛸 수 있습니다.

```bash
python main.py --skip-modeling
```

파이프라인 결과는 아래 경로들에 저장됩니다.

- `data/interim`
- `data/processed`
- `results/plots`
- `results/tables`
- `results/models`
- `results/summaries`

## 참고

- 1차 분석 코호트: `erbb2_status == N` 이고 `pathologic_response_pcr_rd ∈ {pCR, RD}`
- 2차 분석 코호트: `erbb2_status == N` 이고 `pathologic_response_rcb_class ∈ {RCB-0/I, RCB-II, RCB-III}`
- ML의 1차 타깃: `pCR vs RD`
- 2차 타깃은 보조 분석용으로 export만 하며, 기본 모델링 타깃으로는 사용하지 않습니다.
- `GSE25066_RAW.tar`는 후속 민감도 분석용으로만 보관하며, 메인 파이프라인에서는 사용하지 않습니다.

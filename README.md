

## Preview

 **2 Phase chính**:
- **Phase 1 (/src/00_raw_to_curated)**: Filter riêng lẻ từng raw dataset, chuẩn hoá tên , lấy SMILES qua API, MEHC curate
- **Phase 2 (/src/01_build_dataset)**: Merge ADReCS + Pneumotox & SIDER -> gán nhãn từ pattern ADReCS -> remove duplicate -> intergrate CTDBase mechanism -> tạo candidate negatives -> validate negative với PubChem API

---
## Raw Data Description

Raw data nằm trong folder `00_raw_data/` -> Input cho Phase 1

| File thô (Raw File) | Nguồn | Mô tả |
|---|---|---|
| `ADReCS_Respiratory_Drugs.csv` | ADReCS | Chứa các ID thuốc và mô tả tác dụng phụ (ADR_TERM) liên quan đến hô hấp. |
| `pneumotox_full_database.csv` | PNEUMOTOX | Dữ liệu do physician-curated mô tả độc tính phổi, chứa Pattern phân cấp (ví dụ: I.a, II.d). |
| `meddra_all_se.tsv`, `drug_names.tsv` | SIDER  | Database về thuốc và tác dụng phụ. Cung cấp raw data để lọc ra 4 ADR hô hấp mục tiêu. |
| `drugcentral.structures.smiles.tsv` | DrugCentral | Dữ liệu SMILES  của nhiều thuốc -> dùng làm Candidate Negatives. |
| `need2find.csv` | CTDBase | Dữ liệu tương tác hóa chất-bệnh. File này được **copy thủ công từ trang web CTDBase**, chứa danh sách các marker có cơ chế liên quan đến 4 ADR. |

---

## PHASE 1: RAW TO CURATED

Lấy data thô tải từ web/public database, filter độc tính liên quan đến đường hô hấp, gọi API để lấy SMILES, sau đó chạy MEHC - curation

Output: `data/02_curated_raw/`

---

### 1. Nhánh ADReCS & Pneumotox (`01_adrecs`)

| Thông tin | Chi tiết |
|---|---|
| **Raw Input** | `ADReCS_Respiratory_Drugs.csv`, `pneumotox_full_database.csv`  |
| **Output** | `02_curated_raw/adrecs_pneumotox_curated.csv`  |

#### Step 1: `01_adrecs_merge.py` - Merge & Normalize

- Merge 2 dataset ADReCS và Pneumotox thành 1 file duy nhất.
  1. Đọc `ADReCS_Respiratory_Drugs.csv` (cột: `DRUG_ID`, `DRUG_NAME`, `ADR_ID`, `ADR_TERM`, `PubChem_ID`) và `pneumotox_full_database.csv` (cột: `Drug_Name`, `Pattern`).
  2. Normalize tên thuốc: chuyển về lowercase, loại bỏ gốc muối thừa như *hydrochloride, sulfate, sodium, potassium, monohydrate, maleate...*  để tăng khả năng khớp tên giữa 2 nguồn.
  3. Group ADReCS theo tên chuẩn hóa: gộp các ADR_TERM thành chuỗi nối bằng dấu `;`.
  4. Chia 3 nhóm: thuốc có ở cả 2 nguồn (Common), thuốc chỉ có ở ADReCS (ADReCS_only), thuốc chỉ có ở Pneumotox (Pneumotox_only).
  5. **Map PubChem_ID** từ ADReCS sang Pneumotox. Nếu thuốc chỉ có ở Pneumotox  và không có PubChem_ID --> gán `"TO_BE_FOUND"` để xử lí ở bước tiếp theo.
  6. Concat tất cả vào 1 DataFrame, thêm cột `Source` đánh dấu nguồn gốc.
- **Output**: `merged_database_v2.csv` + `adrecs_only_drugs.csv`

#### Step 2: `02_adrecs_fetch_pubchem.py` - Fetch PubChem SMILES

- Lấy  SMILES từ PubChem cho toàn bộ thuốc.
  1.  Gọi PubChem API `/compound/name/{drug_name}/cids/JSON` để tìm PubChem CID từ tên thuốc chưa có CID (có TO_BE_FOUND).
  2. Batch fetch SMILES cho các thuốc đã có CID: Gửi request batch (5 CIDs/request) lên PubChem API `/compound/cid/{cids}/property/CanonicalSMILES/JSON`. Có retry logic (6 lần, backoff exponential). Nếu batch bị lỗi -> fallback sang fetch từng CID một. .
- **Output**: `adrecs_pubchem_smiles.csv`

#### Step 3: `03_adrecs_fetch_chembl.py` - Fetch ChEMBL SMILES (Fallback)

- Bổ sung SMILES cho những thuốc mà PubChem không trả về được. ( 1 số thuốc không có Pubchem ID, nhưng có DrugBank_ID,..)
  1. Lọc ra các dòng có `SMILES` bị rỗng hoặc NaN.
  2. Nếu thuốc có `DrugBank_ID` (bắt đầu bằng "DB") --> gọi ChEMBL API `/chembl/api/data/molecule?molecule_dictionary__drugbank_id={db_id}` để lấy `canonical_smiles`.
- **Output** `adrecs_chembl_smiles.csv`

#### Step 4: `04_adrecs_run_mehc.py` - MEHC Curation
- **Output**: `mehc/refinement/post_refined_smiles.csv`

#### Step 5: `05_adrecs_merge_metadata.py` - Merge Metadata

- Ghép lại metadata gốc (Drug_Name, Pattern, PubChem_ID, Source...) với SMILES đã chuẩn hóa bởi MEHC.
  1. Đọc file raw (trước MEHC) và file curated (sau MEHC).
  2. Drop duplicate theo Drug_Name (giữ bản ghi đầu tiên).
  3. Inner join trên cột `Drug_Name`, chỉ giữ lại những thuốc mà MEHC chấp nhận (loại bỏ thuốc bị MEHC từ chối).
  4. Thay thế cột SMILES cũ bằng SMILES mới từ MEHC.
- **Output cuối**: `02_curated_raw/adrecs_pneumotox_curated.csv` 

---

### 2. Nhánh SIDER (`02_sider`)

| Thông tin | Chi tiết |
|---|---|
| **Raw Input** | `meddra_all_se.tsv`, `drug_names.tsv`  |
| **Output** | `02_curated_raw/sider_respiratory_adrs.csv` |

#### Step 1: `01_sider_fetch_pubchem.py` - Filter Respiratory & Fetch PubChem

- Lọc dataset SIDER để chỉ giữ lại các records liên quan đến 4 bệnh lý hô hấp mục tiêu, đồng thời lấy SMILES.
  1. Đọc `drug_names.tsv`: Tạo mapping `STITCH_ID -> Drug_Name`.
  2. Lọc `meddra_all_se.tsv`: Đọc từng dòng (TSV, cột thứ 6 là side-effect term). Đối chiếu với 4 TARGET_PATTERNS (không phân biệt hoa/thường, loại bỏ khoảng trắng và ký tự đặc biệt).
  3. Gom nhóm theo STITCH_ID: mỗi thuốc sẽ có 1 tập hợp patterns (ví dụ: `"Pleural effusion; Pulmonary oedema"`).
  4. **Chuyển STITCH_ID sang PubChem CID**: Cắt bỏ tiền tố "CID1" để lấy số CID, gọi PubChem API `/compound/cid/{cid}/property/CanonicalSMILES/TXT`. 
- **Output**: `sider_filtered_nomehc.csv`

#### Step 2: `03_sider_prepare_mehc.py` - Format for MEHC

- Prep file cho MEHC - curation
  1. Đọc `sider_filtered_nomehc.csv`.
  2. Chỉ giữ lại 2 cột: `SMILES` và `index` (1-based, làm khóa để ghép metadata về sau).
  3. Loại bỏ toàn bộ metadata (STITCH_ID, Drug_Name, Patterns)
- **Output**: `respiratory_labeled_for_curation.csv`

#### Step 3: `04_sider_run_mehc.py` - MEHC Curation
- **Output**: `mehc/refinement/post_refined_smiles.csv`

#### Step 4: `05_sider_merge_metadata.py` - Merge Metadata
Ghép trả lại metadata gốc cho file SMILES đã curate.
  1. Đọc MEHC output (có cột `smiles`, `index`).
  2. Đọc file raw (có cột `STITCH_ID`, `Drug_Name`, `SMILES`, `Respiratory_Patterns`). Tạo cột `index` 1-based.
  3. Left join trên cột `index` để map metadata về đúng dòng.
  4. Xuất file cuối chỉ gồm 4 cột: `STITCH_ID`, `Drug_Name`, `SMILES`, `Respiratory_Patterns`. 
- **Output**: `02_curated_raw/sider_respiratory_adrs.csv`

---

### 3. Nhánh DrugCentral (`03_drugcentral`)

| Thông tin | Chi tiết |
|---|---|
| **Raw Input** | `drugcentral.structures.smiles.tsv` |
| **Output** | `02_curated_raw/drugcentral_smiles.csv` |

#### Step 1: `01_drugcentral_filter.py` - Filter Safe Drugs

- Filter kho thuốc từ DrugCentral để làm nguồn Negative candidates.
  1. Đọc `pneumotox_adrecs_curated.csv` (kết quả của nhánh Pneumotox + ADReCS) để lấy danh sách tên thuốc đã biết có độc tính.
  2. Đọc `drugcentral.structures.smiles.tsv` (cột: `SMILES`, `InChI`, `InChIKey`, `ID`, `INN`, `CAS_RN`).
  3. **Loại trừ**: So khớp tên thuốc `INN` (lowercase, strip) với danh sách thuốc độc tính từ ADReCS + pneumotox. Giữ lại những thuốc không có trong danh sách độc tính.


- **Output**: `drugcentralonly.csv`

#### Step 2: `02_drugcentral_run_mehc.py` - MEHC Curation
Sau khi chuẩn hóa. Export thẳng ra `02_curated_raw/drugcentral_smiles.csv`.
- **Output**: `02_curated_raw/drugcentral_smiles.csv` 

---

### 4. Nhánh CTDBase (`04_ctdbase`)

| Thông tin | Chi tiết |
|---|---|
| **Raw Input** | `need2find.csv` ( file lọc thủ công từ web) |
| **Output** | `02_curated_raw/ctd_mechanism_markers.csv` |

#### Step 1: `01_ctd_fetch_api.py` - Fetch SMILES
  1. Đọc `need2find.csv` (cột: `Drug_Name`, `Pattern`).
  2. Với mỗi `Drug_Name`, tạo nhiều variant để tăng khả năng khớp:
     - Tên gốc
     - Tách theo dấu phẩy
     - Làm sạch dấu ngoặc (vd: loại bỏ `(...)`)
  3. Thử lần lượt từng variant gọi PubChem API `/compound/name/{name}/property/CanonicalSMILES/TXT`. Trả về SMILES đầu tiên thành công.
- **Output**: `ctdbase_smiles.csv`

#### Step 2: `02_ctd_run_mehc.py` - MEHC Curation

-Gọi MEHC. Export vào thư mục `mehc/refinement/`.
- **Output**: `mehc/refinement/post_refined_smiles.csv`

#### Step 3: `03_ctd_merge_metadata.py` - Merge Metadata

- Ghép cột `Pattern` từ file raw gốc trở lại với SMILES đã chuẩn hóa.
  1. Đọc file raw `need2find.csv` (cột `Drug_Name`, `Pattern`).
  2. Đọc file MEHC output.
  3. **Inner join** trên `Drug_Name` - chỉ giữ những thuốc mà MEHC chấp nhận. - 
  4. output 3 cột: `SMILES`, `Drug_Name`, `Pattern`.
- **Output cuối**: `02_curated_raw/ctd_mechanism_markers.csv`

---

## PHASE 2: BUILD DATASET

Merge 4 file curated, phân tích chuỗi văn bản để gán label, build + validate Negative Pool -> final dataset

Các intermediate được lưu tại: `data/01_intermediate/phase2_assembly/`

---

### Step 1: Merge adrecs + pneumotox + sider -- `01_merge_adrecs_sider.py`

| Script | Action | Output |
|---|---|---|
| `01_merge_adrecs_sider.py` | Concatenate ADReCS/PNEUMOTOX dataset và SIDER dataset | `01_merged_adrecs_sider.csv` |

| Input | Output |
|---|---|
| `adrecs_pneumotox_curated.csv` , `sider_respiratory_adrs.csv`| `01_merged_adrecs_sider.csv`  |

-  Gộp 3 nguồn dữ liệu độc tính phổi chính (ADReCS+Pneumotox và SIDER) thành một pool chung.
  1. Đọc file ADReCS curated (giữ toàn bộ cột).
  2. Đọc file SIDER curated, đổi tên `Respiratory_Patterns` -> `Pattern` để thống nhất tên cột.
  3. `pd.concat()` hai DataFrame, giữ đúng cột của ADReCS làm cơ sở (SIDER thiếu các cột như `PubChem_ID`, `DrugBank_ID` -> tự động điền `NaN`).

---

### Step 2: Label - `02_annotate_labels_from_pattern.py`

| Script | Action | Output |
|---|---|---|
| `02_annotate_labels_from_pattern.py` | Extracts canonical codes (e.g., I.b, V.a, VI.a) from the Pattern column using regex. Also extracts `ADReCS:` text clauses. | `02_labeled_by_pattern.csv` |

| Input | Output |
|---|---|
| `01_merged_adrecs_sider.csv` (1,790) | `02_labeled_by_pattern.csv` (1,790, thêm 4 cột label) |

- Chuyển cột `Pattern` (văn bản mô tả) thành 4 nhãn nhị phân  tương ứng với 4 loại adr

  Mapping rules hơi dài, ở trong script ạ

---

### Step 3: Deduplication by Drug Name -- `03_deduplicate_by_drugname.py`

| Script | Action | Output |
|---|---|---|
| `03_deduplicate_by_drugname.py` | Groups rows by Drug_Name -> resolve conflicts from overlapping sources | `03_deduplicated_by_drugname.csv` |

| Input | Output |
|---|---|
| `02_labeled_by_pattern.csv` | `03_deduplicated_by_drugname.csv` |

- Xóa bỏ các dòng trùng lặp thuốc (cùng Drug_Name nhưng xuất hiện từ nhiều nguồn).
- **Deduplication Logic**: Nếu bản ghi xuất hiện cuối cùng có chứa "ADReCS" trong cột Pattern, script sẽ **giữ lại bản ghi đầu tiên** (thường chứa mô tả chi tiết hơn từ PNEUMOTOX). N
gược lại, bản ghi cuối cùng sẽ được giữ lại. Việc ưu tiên này dựa trên thực tế rằng PNEUMOTOX là cơ sở dữ liệu do bác sĩ biên soạn (physician-curated), trong khi **SIDER** (phát triển tại European Molecular Biology Laboratory - EMBL) là cơ sở dữ liệu toàn diện về thuốc đã lưu hành và ADR được ghi nhận, khiến cả hai nguồn đều đáng tin cậy hơn ADReCS cho việc giải quyết trùng lặp.
 
---

### Step 4: Augmentation with Unlabeled Structures -- `04_augment_drugcentral_smiles.py`

| Script | Action | Output |
|---|---|---|
| `04_augment_drugcentral_smiles.py` | Loads `drugcentral_smiles.csv` and constructs new rows with all 4 labels set to 0 | `04_augmented_drugcentral.csv` |

| Input | Output |
|---|---|
| `03_deduplicated_by_drugname.csv`, `drugcentral_smiles.csv`  | `04_augmented_drugcentral.csv` |

- Nhồi toàn bộ dữ liệu DrugCentral (thuốc an toàn) vào làm **Candidate Negative Pool** (nhãn = 0 cho tất cả 4 cột label).Trước khi add thuốc từ DrugCentral vào làm Negative, đối chiếu và lọc bỏ tất cả những thuốc nào mà Tên của chúng ĐÃ TỒN TẠI trong tập Positive. 
- Bổ sung các phân tử không có bằng chứng trực tiếp về 4 loại ADR mục tiêu.
  1. Đọc file dedup (Positive), đọc file DrugCentral curated.
  2. Tạo DataFrame mới với `Drug_Name = INN`, `Pattern = ''`, tất cả `Label_* = 0`.
  3. `pd.concat()` Positive + DrugCentral candidates.

---

### Step 5: `05_integrate_ctd_markers.py`

| Script | Action | Output |
|---|---|---|
| `05_integrate_ctd_markers.py` | Loads `ctd_mechanism_markers.csv`. The CTD file uses a simplified Pattern column (e.g., plain text "ILD") | `05_merged_ctd_positive_pool.csv` |

| Input | Output |
|---|---|
| `04_augmented_drugcentral.csv` , `ctd_mechanism_markers.csv` | `05_merged_ctd_positive_pool.csv` |

- Đưa CTDBase vào, gán nhãn Positive cho những chất có pattern độc tính rõ ràng.
- **Specific Mapping**:
  - `Pattern == "ILD"` → `Label_ILD = 1`
  - `Pattern == "pulmonary edema"` → `Label_Edema = 1`
  - `Pattern == "pleural effusion"` → `Label_Pleural = 1`
  - `Pattern == "pulmonary embolism"` → `Label_Embolism = 1`
- Bước này bổ sung dữ liệu làm giàu nhóm Positive bằng các marker cơ chế đã được xác nhận sinh học. ( Trên web tìm kiếm theo phenotype thì các thuốc trong list này có liên quan)


---

### Step 6: Validation làm sạch tập negative

#### 6a: Isolation of Candidate Negatives `06a_create_candidate_negatives.py`

| Script | Action | Output |
|---|---|---|
| `06a_create_candidate_negatives.py` | Extracts rows from the Positive Pool where all 4 labels are exactly 0 | `06_candidate_negatives_raw.csv` |

| Input | Output |
|---|---|
| `05_merged_ctd_positive_pool.csv`  | `06_candidate_negatives_raw.csv` |

- Tách toàn bộ các chất có **tổng 4 nhãn = 0** ra thành Candidate Negative Pool riêng để kiểm duyệt.
- **Logic**: `df[df[LABEL_COLS].sum(axis=1) == 0]`

#### 6b: PubChem API Scoring `06b_score_candidates_cautious_cache.py`

| Script | Action | Output |
|---|---|---|
| `06b_score_candidates_cautious_cache.py` | Queries the PubChem PUG-View (for textual summaries) and Assay Summary APIs. Scans only adverse effects sections for curated keywords | `07_pubchem_scored_negatives.csv` |

| Input | Output |
|---|---|
| `06_candidate_negatives_raw.csv` | `07_pubchem_scored_negatives.csv` |

-  validate từng Candidate Negative bằng PubChem để đảm bảo không có chất nào thực sự có độc tính phổi bị lọt vào nhóm Negative
- **Scoring Details**:
  1. **Resolve CID**: Thử 3 chiều -- `PubChem_ID` -> `Drug_Name`  -> `Drug_Name` 
  2. **Chấm điểm qua PUG-View**: Gọi PubChem PUG-View API để lấy toàn bộ thông tin về hợp chất. Quét chỉ các section có tiêu đề chứa **"Adverse Effects"**, **"Toxicity"**, **"Safety"**, và **"Warning"** cho danh sách keyword đã compile cho 4 bệnh mục tiêu (30+ từ khóa bao gồm: `interstitial lung disease, pulmonary fibrosis, pneumonitis, pleural effusion, pulmonary embolism, deep vein thrombosis`...).
  3. **Chấm điểm qua Assay Summary**: Tương tự, gọi Assay Summary API, chỉ đếm keyword trong các assay có chứa adverse terms.
  4. **Total Score = (Keyword matches in PUG-View) + (Keyword matches in Assay Summary)**. Nếu `score > 0` --> thuốc này có dấu hiệu độc tính phổi --> **loại bỏ** khỏi Negative Pool.


> Warning: Chạy lâu 

> Note: "" "Assay Summary là dịch vụ của PubChem cung cấp bản tóm tắt mô tả của một BioAssay (thí nghiệm sinh học), được truy cập qua AID (Assay ID). PUG-View (Power User Gateway View) là một dịch vụ web REST-style của PubChem, chuyên dùng để truy cập dữ liệu chú thích (annotation) được tích hợp trong các bản ghi PubChem – tức là những thông tin đã được tổng hợp, chuẩn hóa từ nhiều nguồn khác nhau (nhà sản xuất, cơ quan quản lý, tài liệu khoa học) """
---

### Step 7: Final Dataset Assembly -- `07_assemble_final_dataset.py`

| Script | Action | Output |
|---|---|---|
| `07_assemble_final_dataset.py` | Constructs the final balanced dataset by concatenating positive and negative sets | `balanced_pulmonary_adr_dataset.csv` |

| Input | Output |
|---|---|
| `05_merged_ctd_positive_pool.csv`, `07_pubchem_scored_negatives.csv` | `balanced_pulmonary_adr_dataset.csv`  |

-  Ghép Positive Pool + Verified Negative Pool thành dataset cuối cùng.
  1. **Positive Set**: Tất cả các dòng từ `05_merged_ctd_positive_pool.csv` có **ít nhất 1 nhãn == 1** 

  2. **Negative Set**: Tất cả các dòng từ `07_pubchem_scored_negatives.csv` có **score == 0** (tức là không có bất kỳ bằng chứng nào về 4 loại ADR trong toàn bộ PubChem database) 
  3. `pd.concat()` Positive + Negative.
  4. **Dedup theo SMILES** (`drop_duplicates(subset=['SMILES'], keep='first')`) để loại bỏ các dòng có cùng cấu trúc hóa học.
  5. Sắp xếp cột: `Drug_Name, Pattern, PubChem_ID, DrugBank_ID, KEGG_ID, SMILES, Label_ILD, Label_Edema, Label_Pleural, Label_Embolism`.


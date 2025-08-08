# WebAR Cloud Model Generator

このプロジェクトは、気象衛星ひまわりのTBB（等価黒体輝度温度）データを読み込み、指定領域の雲の3Dモデル（glb形式）を生成するPythonスクリプトです。AR.js + A-Frameを使ったマーカーベースのWebARでの表示を目的としており、HTMLファイルのサンプルも提供します。

## ■ サンプル

- [AR技術を活用した立体感覚の涵養を目的とする地学教材の開発／雲の立体モデル（AR教材のサンプル）](https://robo.mydns.jp/WebAR/index.html#kumonoRittai)

## ■ 必要な環境

- Python 3.8 以上
- 利用するモジュール:
  `numpy`, `xarray`, `trimesh`, `geopandas`, `shapely`
- `trimesh`の依存ライブラリ（`scipy`など）も必要です。

## ■ インストール手順

```bash
git clone https://github.com/your-username/CloudModelGenerator.git
cd CloudModelGenerator
pip install numpy xarray trimesh geopandas shapely
```

## ■ 使用方法

### 1. 海岸線データの準備

このスクリプトは、雲モデルと共に表示する海岸線データを必要とします。[Natural Earth](https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/) から10m解像度の海岸線データ（`ne_10m_coastline.zip`）をダウンロードし、解凍して得られる `ne_10m_coastline.geojson` をスクリプトと同じディレクトリに配置してください。

### 2. 雲データの取得とGLBモデルの生成

TBBデータは、千葉大学環境リモートセンシング研究センターのFTPサーバーから自動的にダウンロードされます。

以下のコマンドを実行して、指定した領域と日時の雲モデルを生成します。

```bash
python3 makeCloudGLB.py REGION YYYYMMDDHH
```

- **REGION**: 表示する領域を指定します。以下のいずれかを選択できます。
  - `japan`: 日本全国
  - `southjapan`: 西日本
  - `hokkaido`: 北海道
  - `kantou`: 関東
- **YYYYMMDDHH**: 世界標準時（UTC）での年月日と時間を指定します。（例: `2024082800`）

実行が成功すると、`CM_REGION_YYYYMMDDHH.glb` という名前のファイルが生成されます。

### 3. WebAR で表示する

生成された `.glb` ファイルと `sample.html` を同じディレクトリに置き、Webサーバーで公開します。その際、`sample.html` 内の`clouds.glb` を生成された `.glb` ファイルのファイル名に置き換えます。また、scale="0.003 0.030 0.003" がモデルを表示する際の倍率を表しています。適宜調整してください。２つめが鉛直方向のスケールです。ファイルの設定では、水平スケールに対して鉛直方向は10倍に拡大されています。

スマートフォンなどで `sample.html` にアクセスし、マーカー（Hiroなど）をカメラに映せば、AR上に雲の立体モデルが表示されます。

## ■ データの出典

- 雲データ: [千葉大学 環境リモートセンシング研究センター](http://www.cr.chiba-u.jp/databases/gridded-data-jp.html)
- 海岸線データ: [Natural Earth](https://www.naturalearthdata.com/)

---

## ■ ファイル構成

```
CloudModelGenerator/
├── makeCloudGLB.py              # メインスクリプト
├── sample.html                  # A-Frame + AR.js による表示例
├── ne_10m_coastline.geojson     # 海岸線データ
└── README.md                    # このファイル
```

---

## 📝 ライセンス

このプロジェクトは MIT ライセンスのもとで公開されています。
